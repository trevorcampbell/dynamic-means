#ifndef __KERNDYNMEANS_IMPL_HPP

template<typename G>
KernDynMeans<G>::KernDynMeans(double lambda, double Q, double tau, bool verbose){
	if (lambda < 0 || Q < 0 || tau < 0){
		cout << "libkerndynmeans: ERROR: Parameters of Kernel Dynamic Means cannot be < 0." << endl;
		cout << "libkerndynmeans: Lambda: " << lambda << " Q: " << Q << " tau: " << tau << endl;
	}
	this->verbose = verbose;
	this->maxLblPrevUsed = -1;
	this->ages.clear();
	this->oldprmlbls.clear();
	this->weights.clear();
	this->gammas.clear();
	this->agecosts.clear();
	this->lambda = lambda;
	this->Q = Q;
	this->tau = tau;
	this->sigma = this->sigmaUB = this->sigmaLB = 0;
	try{
		this->grbenv = new GRBEnv();
		//start up GRB
		grbenv->set(GRB_IntParam_OutputFlag, 0);//controls the output of Gurobi - 0 means no output, 1 means normal output
		//grbenv->set(GRB_IntParam_Method, 1);//controls which method Gurobi uses - default -1 (auto), 0=primal simplex, 1=dual simplex, 2=barrier, 3=concurrent, 4=deterministic concurrent
		grbenv->set(GRB_IntParam_Threads, 1);//controls the number of threads Gurobi uses - I force it to use 1 since the optimization in this algorithm is fairly small/simple
										     	//											and it just ends up wasting time constantly creating/deleting threads otherwise
	} catch (GRBException e){
		cout << e.getErrorCode() << " " << e.getMessage() << endl;
	}
}

template<typename G>
KernDynMeans<G>::~KernDynMeans(){
	delete this->grbenv;
}

template<typename G>
void KernDynMeans<G>::reset(){
	this->maxLblPrevUsed = -1;
	this->ages.clear();
	this->oldprmlbls.clear();
	this->weights.clear();
	this->gammas.clear();
	this->agecosts.clear();
	this->sigma = this->sigmaUB = this->sigmaLB = 0;
}

//This function updates the weights/ages of all the clusters after each clustering step is complete
template <typename G>
void KernDynMeans<G>::finalizeStep(const G& aff, const vector<int>& lbls, vector<double>& prevgammas_out, vector<int>& prmlbls_out){
	//first increment the age of everything
	//this will be undone below for any current clusters
	for (int i = 0; i < this->ages.size(); i++){
		this->ages[i]++;
	}
	//get the set of current cluster labels
	vector<int> unqlbls = lbls;
	sort(unqlbls.begin(), unqlbls.end());
	unqlbls.erase(unique(unqlbls.begin(), unqlbls.end()), unqlbls.end());
	//for every current cluster
	for (int i = 0; i < unqlbls.size(); i++){
		const int& lbl = unqlbls[i];
		//get the number of nodes assigned to the cluster
		int nclus = 0;
		for (int j = 0; j < lbls.size(); j++){
			nclus += (lbls[j] == lbl ? aff.getNodeCt(j) : 0);
		}
		//check if it's an old label
		auto it2 = find(this->oldprmlbls.begin(), this->oldprmlbls.end(), lbl);
		if (it2 == this->oldprmlbls.end()){
			//it's a new cluster, so create stuff for it
			this->ages.push_back(1);
			this->weights.push_back(nclus);
			this->oldprmlbls.push_back(lbl);
		} else {
			//it's an old cluster, so update old stuff
			int oldidx = distance(this->oldprmlbls.begin(), it2);
			this->ages[oldidx] = 1;
			this->weights[oldidx] = this->gammas[oldidx] + nclus;
		}
	}
	//save the old gammas for output, update this->gammas
	//prevgammas_out still needs to be cleaned up and modified a bit (occurs in the rest of this function)
	prevgammas_out = this->gammas;
	//pad prevgammas with zeros for all new clusters
	prevgammas_out.insert(prevgammas_out.end(), this->ages.size()-prevgammas_out.size(), 0);
	//update the internal gammas
	this->gammas.clear();
	this->gammas.reserve(this->ages.size());
	for (int i = 0; i < this->ages.size(); i++){
			this->gammas.push_back(1.0/(1.0/this->weights[i] + this->ages[i]*this->tau)); 
	}
	//update the Age Costs
	this->agecosts.clear();
	this->agecosts.reserve(this->ages.size());
	for (int i = 0; i < this->ages.size(); i++){
			this->agecosts.push_back(this->Q*this->ages[i]); 
	}

	//delete any old cluster whose age cost exceeds lambda
	for (int i = 0; i < this->ages.size(); i++){
		if (this->agecosts[i] > this->lambda){
			this->weights.erase(this->weights.begin()+i);
			this->ages.erase(this->ages.begin()+i);
			this->gammas.erase(this->gammas.begin()+i);
			prevgammas_out.erase(prevgammas_out.begin()+i);
			this->agecosts.erase(this->agecosts.begin()+i);
			this->oldprmlbls.erase(this->oldprmlbls.begin()+i);
			i--;
		}
	}
	prmlbls_out = this->oldprmlbls;//save the parameter labels output 

	//update this->maxLblPrevUsed if new clusters were created
	int maxlbl = *max_element(lbls.begin(), lbls.end());
	if (maxlbl >= this->maxLblPrevUsed){
		this->maxLblPrevUsed = maxlbl;
	}

	return;
}

template<typename G>
void KernDynMeans<G>::cluster(const G& aff, const int nRestarts, const int nCoarsest, std::vector<int>& finalLabels, double& finalObj, std::vector<double>& finalGammas, 
		std::vector<int>& finalPrmLbls, double& tTaken){
	timeval tStart;
	gettimeofday(&tStart, NULL);

	int nNodes = aff.getNNodes();
	int nOldPrms = this->oldprmlbls.size();

	if (nNodes <= 0){
		cout << "libkerndynmeans: WARNING: nNodes <=0 (= " << nNodes << "); Returning empty labels."<<  endl;
		finalLabels = vector<int>();
		tTaken = 0;
		finalObj = 0;
		return;
	}
	if (nRestarts <= 0){
		cout << "libkerndynmeans: ERROR: nRestarts <=0 (= " << nRestarts << ")"<<  endl;
		return;
	}
	//compute sigma upper bound via diagonal dominance and initialize sigma to 0
	if (verbose){
		cout << "libkerndynmeans: Computing sigma bounds." << endl;
	}
	this->initializeSigma(aff);
	if (verbose){
		cout << "libkerndynmeans: Sigma bounds: [" << this->sigmaLB << ", " << this->sigmaUB << "]." << endl;
		cout << "libkerndynmeans: Clustering " << nNodes << " datapoints with " << nRestarts << " restarts." << endl;
		cout << "libkerndynmeans: " << nOldPrms << " old clusters possibly alive from previous timesteps." << endl;
	}

	std::vector<int> minLbls;
	double minObj = std::numeric_limits<double>::max();
	for(int rest = 0; rest < nRestarts; rest++){
		if (verbose){
			cout << "libkerndynmeans: Attempt " << rest+1 << "/" << nRestarts << ", Minimum Obj = " << minObj << endl;
		}

		//first, form the coarsification levels in the graph if necessary
		std::vector<int> lbls;
		if(nNodes > nCoarsest){
			if (verbose){
				cout << "libkerndynmeans: Coarsifying " << nNodes << " nodes at level 0." << endl;
			}
			std::stack<CoarseGraph<G> > coarsestack; //stores coarsified graphs
			CoarseGraph<G> cg;
			cg.coarsify(aff);
			coarsestack.push(cg);
			while(coarsestack.top().getNNodes() > nCoarsest){
				if (verbose){
					cout << "libkerndynmeans: Coarsifying " << coarsestack.top().getNNodes() << " nodes at level " << coarsestack.size() << "." << endl;
				}
				CoarseGraph<G> cg2;
				cg2.coarsify(coarsestack.top());
				coarsestack.push(cg2);
			}
			if (verbose){
				cout << "libkerndynmeans: Done coarsifying, top level " << coarsestack.size() << " has " << coarsestack.top().getNNodes() << " nodes." << endl;
			}
			//next, step down through the refinements and cluster, initializing from the coarser level
			while(!coarsestack.empty()){
				if (verbose){
					cout << "libkerndynmeans: Running clustering at level " << coarsestack.size() << " with " << coarsestack.top().getNNodes() << " nodes." << endl;
				}
				//optimize the labels for the current top of coarsestack
				lbls = this->clusterAtLevel(coarsestack.top(), lbls); //lbls starts out empty, clusterAtLevel knows to use a base clustering
				//refine the labels
				lbls = coarsestack.top().getRefinedLabels(lbls);
				coarsestack.pop();
			}
		}
		if (verbose){
			cout << "libkerndynmeans: Running clustering at data level." << endl;
		}
		//final clustering at the data level
		lbls = this->clusterAtLevel(aff, lbls);

		//finally, compute the kernelized dynamic means objective
		if (verbose){
			cout << "libkerndynmeans: Computing objective..." << endl;
		}
		double obj = this->objective(aff, lbls);
		if (verbose){
			cout << "libkerndynmeans: Objective = " << obj << endl;
		}
		if (obj < minObj){
			minLbls = lbls;
			minObj = obj;
		}
	}

	if (verbose){
		vector<int> unqlbls = minLbls;
		sort(unqlbls.begin(), unqlbls.end());
		unqlbls.erase(unique(unqlbls.begin(), unqlbls.end()), unqlbls.end());
		int numnew = 0;
		for (int i = 0; i < unqlbls.size(); i++){
			if (find(this->oldprmlbls.begin(), this->oldprmlbls.end(), unqlbls[i]) == this->oldprmlbls.end()){
				numnew++;
			}
		}
		int numoldinst = unqlbls.size() - numnew;
		int numolduninst = this->ages.size() - numoldinst;
		cout << endl << "libkerndynmeans: Done clustering. Min Objective: " << minObj << " Old Uninst: " << numolduninst  << " Old Inst: " << numoldinst  << " New: " << numnew <<  endl;
	}
	//update the state of the ddp chain
	this->finalizeStep(aff, minLbls, finalGammas, finalPrmLbls);
	//collect results
	finalObj =  minObj;
	finalLabels = minLbls;
	timeval tCur;
	gettimeofday(&tCur, NULL);
	tTaken = (double)(tCur.tv_sec - tStart.tv_sec) + (double)(tCur.tv_usec - tStart.tv_usec)/1.0e6;
	return;
}


template<typename G>
template <typename T> 
std::vector<int> KernDynMeans<G>::clusterAtLevel(const T& aff, std::vector<int> lbls){ 
	if (lbls.size() < aff.getNNodes()){ // Base Clustering -- Use spectral clustering on data, maximum bipartite matching to link old clusters
		if (verbose){ cout << "Running base spectral clustering..." << endl;}
		//get the data labels from spectral clustering
		lbls = this->baseCluster(aff);
		//find the optimal correspondence between old/current clusters
		lbls = this->updateOldNewCorrespondence(aff, lbls);
		//initlbls is now ready for regular refinement iterations
		if (verbose){ cout << "Done base spectral clustering with objective: " << this->objective(aff, lbls) << endl;}
	}

	//run the refinement iterations
	double prevobj = this->objective(aff, lbls);
	double diff = 1.0;
	int itr = 0;
	while(diff > 1e-6){
		this->sigma = this->sigmaLB; //start sigma at its lower bound
		itr++;
		std::vector<int> tmplbls = this->updateLabels(aff, lbls);
		double tmpobj = this->objective(aff, tmplbls);
		//if the update increased the objective, update the sigma lower bound by searching backwards from sigmaub
		if (tmpobj > prevobj && fabs(tmpobj-prevobj)>1e-6){ //1e-6 to help get rid of numerical noise
			if (verbose){
				cout << "libkerndynmeans: Monotonicity violated!" << endl;
				cout << "libkerndynmeans: Finding new sigmaLB..." << endl;
			}
			this->sigma = this->sigmaUB; //start at the upper bound
			do{
				this->sigma = (this->sigma + this->sigmaLB)/2.0; //progressively backwards search towards LB
				tmplbls = this->updateLabels(aff, lbls);
				tmpobj = this->objective(aff, tmplbls);
			} while (tmpobj < prevobj); //if we find the sigma that violates monotonicity
			this->sigmaLB += 2.0*(this->sigma-this->sigmaLB); //set sigmaLB to the last one that worked -- this is a conservative lower bound so won't cause cycling
			if (verbose){
				cout << "libkerndynmeans: New sigma bounds are [" << this->sigmaLB << ", " << this->sigmaUB << "]." << endl;
			}
		}
		lbls = this->updateOldNewCorrespondence(aff, tmplbls); //guaranteed not to increase objective, no check needed
		double obj = this->objective(aff, lbls);
		diff = fabs((obj-prevobj)/obj);
		prevobj = obj;
		if (verbose){ cout << "libkerndynmeans: Kernelized clustering iteration " << itr << ", obj = " << obj << endl;}
	}
	if (verbose){cout << endl;}
	return lbls;
}


template <typename G>
template <typename T>
std::vector<int> KernDynMeans<G>::updateLabels(const T& aff, std::vector<int> lbls) const{

	//get the unique labels
	vector<int> unqlbls = lbls;
	sort(unqlbls.begin(), unqlbls.end());
	unqlbls.erase(unique(unqlbls.begin(), unqlbls.end()), unqlbls.end());


	//get the observations in each cluster
	//and the sizes of each cluster
	map<int, std::vector<int> > idcsInClus;
	map<int, double > nInClus;
	for (int i = 0; i < lbls.size(); i++){
		idcsInClus[lbls[i]].push_back(i);
		if(nInClus.count(lbls[i]) == 0){
			nInClus[lbls[i]] = 0.0;
		}
		nInClus[lbls[i]] += aff.getNodeCt(i);
	}

	//precompute the within-cluster sums
	std::map<int, double> inClusterSum;
	for (int i = 0; i < unqlbls.size(); i++){
		double sum = 0;
		const std::vector<int>& clus = idcsInClus[unqlbls[i]];
		const int& lbl = unqlbls[i];
		for (int k = 0; k < clus.size(); k++){
			int nct = aff.getNodeCt(clus[k]);
			sum += aff.diagSelfSimDD(clus[k])+nct*this->sigma + 2.0*aff.offDiagSelfSimDD(clus[k]);
			for (int m = k+1; m < clus.size(); m++){
				sum += 2.0*aff.simDD(clus[k], clus[m]);
			}
		}
		inClusterSum[lbl] = sum;
	}

	//minimize the cost associated with each observation individually based on the old labelling
	std::vector<int> newlbls = lbls;
	int nextlbl = std::max(this->maxLblPrevUsed+1, 1+ *max_element(lbls.begin(), lbls.end()));//for this round, handles labelling of new clusters
	//need max_element and maxlblprevused above to avoid assigning new clustered data to old uninstantiated clusters or current new clusters from the previous iteration of labellign
	for (int i = 0; i < lbls.size(); i++){
		int nct = aff.getNodeCt(i);
		double minCost = this->lambda + this->sigma + (1.0-1.0/nct)*(aff.diagSelfSimDD(i) + nct*this->sigma) 
							- 2.0/nct*aff.offDiagSelfSimDD(i); //default to creating a new cluster, and then try to beat it 
		int minLbl = -1;
		const int& prevlbl = lbls[i];

		//if there's only one node in the cluster, and no earlier observation was assigned to it, remove it before proceeding
		if (idcsInClus[prevlbl].size() == 1 && find(newlbls.begin(), newlbls.begin()+i, prevlbl) == (newlbls.begin()+i)){
			unqlbls.erase(std::remove(unqlbls.begin(), unqlbls.end(), prevlbl), unqlbls.end());
			idcsInClus.erase(prevlbl);
			nInClus.erase(prevlbl);
			inClusterSum.erase(prevlbl);
		}

		//run through instantiated clusters
		for (int k = 0; k < unqlbls.size(); k++){
			const std::vector<int>& clus = idcsInClus[unqlbls[k]];
			const int& lbl = unqlbls[k];
			double cost = 0;
			auto it = find(this->oldprmlbls.begin(), this->oldprmlbls.end(), lbl);
			if (it == this->oldprmlbls.end()){//if it's a new cluster in this timestep, no gamma stuff is needed
				cost = aff.diagSelfSimDD(i)+nct*this->sigma+(double)nct/(nInClus[lbl]*nInClus[lbl])*inClusterSum[lbl];
				if (prevlbl == lbl){
					cost += -2.0/nInClus[lbl]*(aff.diagSelfSimDD(i)+nct*this->sigma+2.0*aff.offDiagSelfSimDD(i));
				}
				for (int j = 0; j < clus.size(); j++){
					if (clus[j] != i){
						cost += -2.0/nInClus[lbl]*aff.simDD(i, clus[j]);
					}
				}
			} else {//it's an instantiated old cluster, need to do gamma stuff
				double oldidx = distance(this->oldprmlbls.begin(), it);
				double factor = 1.0/(nInClus[lbl]+this->gammas[oldidx]);
				cost = aff.diagSelfSimDD(i) +nct*this->sigma
					-2.0*this->gammas[oldidx]*factor*aff.simDP(i, oldidx) 
					+ (double)nct*this->gammas[oldidx]*this->gammas[oldidx]*factor*factor*(aff.selfSimPP(oldidx) + this->sigma/this->gammas[oldidx])
					+(double)nct*factor*factor*inClusterSum[lbl];
				if (prevlbl == lbl){
					//need to be careful about similarities when the observation was previously in this cluster
					cost += -2.0*factor*(aff.diagSelfSimDD(i)+nct*this->sigma+2.0*aff.offDiagSelfSimDD(i));
				}
				for (int j = 0; j < clus.size(); j++){
					cost += 2.0*factor*factor*nct*this->gammas[oldidx]*aff.simDP(clus[j], oldidx);
					if (clus[j] != i){
						cost += -2.0*factor*aff.simDD(i, clus[j]);
					}
				}
			}
			if (cost < minCost){
				minCost = cost;
				minLbl = lbl;
			}
		}
		//run through old uninstantiated clusters
		for (int k = 0; k < this->oldprmlbls.size(); k++){
			auto it = find(unqlbls.begin(), unqlbls.end(), this->oldprmlbls[k]);
			if (it == unqlbls.end()){
				double cost = this->agecosts[k]
						+(1.0-1.0/(this->gammas[k]+nct))*(aff.diagSelfSimDD(i)+nct*this->sigma)
						+this->gammas[k]*nct/(this->gammas[k]+nct)*(aff.selfSimPP(k) + this->sigma/this->gammas[k])
						-2.0/(this->gammas[k]+nct)*aff.offDiagSelfSimDD(i)
						-2.0*this->gammas[k]/(this->gammas[k]+nct)*aff.simDP(i, k);
				if (cost < minCost){
					minCost = cost;
					minLbl = this->oldprmlbls[k];
				}
			}
		}
		if (minLbl == -1){
			//create a new cluster
			minLbl = nextlbl;
			nextlbl++;
		}
		//relabel the datapoint
		newlbls[i] = minLbl;
		//if a new cluster was created or old one revived, add stuff into unq/idcsinclus/inclustersum/nInClus/etc
		if (find(unqlbls.begin(), unqlbls.end(), minLbl) == unqlbls.end()){
			unqlbls.push_back(minLbl);
			idcsInClus[minLbl].push_back(i);   //no problem with duplicating the datapoint here (it will still exist in idcsInClus in the old cluster)
												//this is because the distance comparisons to previous cluster centers shouldn't be affected by creating new centers
												//ergo, leave the old dInClus alone, but create a new entry so that other observations can switch tothis cluster
			inClusterSum[minLbl] = aff.diagSelfSimDD(i)+nct*this->sigma+2.0*aff.offDiagSelfSimDD(i);
			nInClus[minLbl] = nct;
		}
	}

	return newlbls;
}


template <typename G>
template <typename T> 
std::vector<int> KernDynMeans<G>::updateOldNewCorrespondence(const T& aff, std::vector<int> lbls) const{
	//get the unique labels
	vector<int> unqlbls = lbls;
	sort(unqlbls.begin(), unqlbls.end());
	unqlbls.erase(unique(unqlbls.begin(), unqlbls.end()), unqlbls.end());

	//get the observations in each cluster
	//and the sizes of each cluster
	map<int, std::vector<int> > idcsInClus;
	map<int, double > nInClus;
	for (int i = 0; i < lbls.size(); i++){
		idcsInClus[lbls[i]].push_back(i);
		if(nInClus.count(lbls[i]) == 0){
			nInClus[lbls[i]] = 0.0;
		}
		nInClus[lbls[i]] += aff.getNodeCt(i);
	}

	//compute the squared cluster sums
	std::map<int, double> inClusterSum;
	for (int i = 0; i < unqlbls.size(); i++){
		double sum = 0;
		const std::vector<int>& clus = idcsInClus[unqlbls[i]];
		const int& lbl = unqlbls[i];
		for (int k = 0; k < clus.size(); k++){
			sum += aff.diagSelfSimDD(clus[k]) + 2.0*aff.offDiagSelfSimDD(clus[k]);
			for (int m = k+1; m < clus.size(); m++){
				sum += 2.0*aff.simDD(clus[k], clus[m]);
			}
		}
		inClusterSum[lbl] = sum;
	}

	//get the old/new correspondences from bipartite matching
 	vector< pair<int, int> > nodePairs; //new clusters in .first, old clusters + one null cluster in .second
 	vector< double > edgeWeights;
	for (int i = 0; i < unqlbls.size(); i++){
		const std::vector<int>& clus = idcsInClus[unqlbls[i]];
		const int& lbl = unqlbls[i];
		for (int j = 0; j < this->oldprmlbls.size(); j++){
			double ewt = this->agecosts[j]
						+ this->gammas[j]*nInClus[lbl]/(this->gammas[j]+nInClus[lbl])*aff.selfSimPP(j)
						-1.0/(this->gammas[j]+nInClus[lbl])*inClusterSum[lbl];
			double oldClusSum = 0;
			for (int k = 0; k < clus.size(); k++){
				oldClusSum += aff.simDP(clus[k], j);
			}
			ewt += -2.0*this->gammas[j]/(this->gammas[j]+nInClus[lbl])*oldClusSum;
			nodePairs.push_back(std::pair<int, int>(lbl, this->oldprmlbls[j]) );
			edgeWeights.push_back(ewt);
		}
		//-1 is the new cluster option
		nodePairs.push_back( std::pair<int, int>(lbl, -1) );
		edgeWeights.push_back(this->lambda-1.0/nInClus[lbl]*inClusterSum[lbl]);
	}
	map<int, int> matching = this->getMinWtMatching(nodePairs, edgeWeights);

	//relabel lbls based on the old/new correspondences
	std::vector<int> newlbls = lbls;
	int nextlbl = this->maxLblPrevUsed+1;
	for (auto it = matching.begin(); it != matching.end(); ++it){
		if (it->second != -1){ //if the current cluster isn't new
			//replace all labels in lbls to the old cluster label
			for (int i = 0; i < lbls.size(); i++){
				if (lbls[i] == it->first){
					newlbls[i] = it->second;
				}
			}
		} else {
			//replace all labels in lbls to the new cluster label
			for (int i = 0; i < lbls.size(); i++){
				if (lbls[i] == it->first){
					newlbls[i] = nextlbl;
				}
			}
			nextlbl++;
		}
	}
	return newlbls;
}

template <typename G>
map<int, int> KernDynMeans<G>::getMinWtMatching(vector< pair<int, int> > nodePairs, vector<double> edgeWeights ) const{
	//used to get the matching between old/current clusters
	//for each  current cluster (on the left, A), there must be one edge outgoing from it
	//for each old cluster (on the right, B), there can be at most one edge outgoing
	//for the final -1 node (on the right, in B), there can be as many outgoing edges as necessary


	//get params
	int nVars = edgeWeights.size();

	try{
	GRBModel grbmodel(*grbenv);
	//add variables/objective
	double* obj = new double[nVars];
	for (int i = 0; i < nVars; i++){
		obj[i] = edgeWeights[i];
	}
	GRBVar* grbvars = grbmodel.addVars(NULL, NULL,obj, NULL, NULL, nVars);

	grbmodel.update();

	//one constraint for each A/B node, plus one constraint for each edge
	vector<int> A, B;
	for (int i = 0; i < nodePairs.size(); i++){
		if (find(A.begin(), A.end(), nodePairs[i].first) == A.end()){
			A.push_back(nodePairs[i].first);
		}
		if (find(B.begin(), B.end(), nodePairs[i].second) == B.end()){
			B.push_back(nodePairs[i].second);
		}
	}
	//add constraints
	//constraint type 1: sum of outgoing edges from A nodes = 1
	for (int i = 0; i < A.size(); i++){
		GRBLinExpr constrlhs;
		for (int j = 0; j < nVars; j++){
			if (nodePairs[j].first == A[i]){
				constrlhs += 1.0*grbvars[j];
			}
		}
		grbmodel.addConstr(constrlhs, GRB_EQUAL, 1);
	}
	//constraint type 2: sum of incoming edges to B nodes <= 1
	//if B[i] == -1, no constraint -- can have as many outgoing edges as required
	for (int i = 0; i < B.size(); i++){
		if (B[i] != -1){
			GRBLinExpr constrlhs;
			for (int j = 0; j < nVars; j++){
				if (nodePairs[j].second == B[i]){
					constrlhs += 1.0*grbvars[j];
				}
			}
			grbmodel.addConstr(constrlhs, GRB_LESS_EQUAL, 1);
		}
	}
	//constraint type 3: all edge variables >= 0
	//the polytope has an implicit bound of >=0 on all variables, don't need this

	//defaults to minimization
	grbmodel.optimize();

	map<int, int> retmap;
	for (int j = 0; j < edgeWeights.size(); j++){
		double val = grbvars[j].get(GRB_DoubleAttr_X);
		if (fabs(val - 1.0) < 1e-10){
			retmap[nodePairs[j].first] = nodePairs[j].second;
			A.erase(remove(A.begin(), A.end(), nodePairs[j].first), A.end());
		}
	}
	delete[] grbvars;
	delete[] obj;
	return retmap;
	} catch (GRBException e){
		cout << "libkerndynmeans: ERROR: Gurobi Error code = " << e.getErrorCode() << endl;
		cout << "libkerndynmeans: MESSAGE: " << e.getMessage() << endl;
	} catch (...){
		cout << "libkerndynmeans: ERROR: Unhandled Gurobi exception during optimization" << endl;
	}
}

template<typename G>
template<typename T> 
double KernDynMeans<G>::objective(const T& aff, const std::vector<int>& lbls) const{
	double cost = 0;
	//get the observations in each cluster
	//and the sizes of each cluster
	map<int, std::vector<int> > idcsInClus;
	map<int, double > nInClus;
	for (int i = 0; i < lbls.size(); i++){
		idcsInClus[lbls[i]].push_back(i);
		if(nInClus.count(lbls[i]) == 0){
			nInClus[lbls[i]] = 0.0;
		}
		nInClus[lbls[i]] += aff.getNodeCt(i);
	}
	//for every current cluster
	for (auto it = idcsInClus.begin(); it != idcsInClus.end(); ++it){
		const int& lbl = it->first;
		const std::vector<int>& clus = it->second;
		//check if it's an old label
		auto it2 = find(this->oldprmlbls.begin(), this->oldprmlbls.end(), lbl);
		if (it2 == this->oldprmlbls.end()){ //it's a new cluster
			cost += this->lambda;//new cluster penalty
			//ratio association term
			for (int i = 0; i < clus.size(); i++){
				cost += (1.0-1.0/nInClus[lbl])*aff.diagSelfSimDD(clus[i]) - 2.0/nInClus[lbl]*aff.offDiagSelfSimDD(clus[i]);
				for (int j = i+1; j < clus.size(); j++){
					cost += -2.0/nInClus[lbl]*aff.simDD(clus[i], clus[j]);
				}
			}
		} else { //it's an old cluster
			int oldidx = distance(this->oldprmlbls.begin(), it2);
			cost += this->agecosts[oldidx];//old cluster penalty
			//ratio association term
			for (int i = 0; i < clus.size(); i++){
				cost += (1.0 - 1.0/(this->gammas[oldidx]+nInClus[lbl]))*aff.diagSelfSimDD(clus[i]) - 2.0/(this->gammas[oldidx]+nInClus[lbl])*aff.offDiagSelfSimDD(clus[i]);
				for (int j = i+1; j < clus.size(); j++){
					cost += -2.0/(this->gammas[oldidx]+nInClus[lbl])*aff.simDD(clus[i], clus[j]);
				}
			}
			cost += this->gammas[oldidx]*nInClus[lbl]/(this->gammas[oldidx]+nInClus[lbl])*aff.selfSimPP(oldidx);
			//old prm ratio association term
			for (int i = 0; i < clus.size(); i++){
				cost += -2.0*this->gammas[oldidx]/(this->gammas[oldidx]+nInClus[lbl])*aff.simDP(clus[i], oldidx);
			}
		}
	}
	return cost;
}

template<typename G>
template<typename T>
std::vector<int> KernDynMeans<G>::baseCluster(const T& aff) const{
	//compute the kernel matrix
	int nA = aff.getNNodes();
	MXd K = MXd::Zero(nA, nA);
	for (int i = 0; i < nA; i++){
		K(i, i) = aff.diagSelfSimDD(i) + 2.0*aff.offDiagSelfSimDD(i);
		for (int j = i+1; j < nA; j++){
			K(i, j) = aff.simDD(i, j);
			K(j, i) = K(i, j);
		}
	}
	//solve the eigensystem for eigenvectors
	Eigen::SelfAdjointEigenSolver<MXd> eigsol;
	eigsol.compute(K);
	//since the eigenvalues are sorted in increasing order, chop off the ones at the front
	VXd eigvals = eigsol.eigenvalues();
	MXd Z = eigsol.eigenvectors();
	int chopIdx = 0;
	while (chopIdx < eigvals.size() && eigvals(chopIdx) < this->lambda) chopIdx++; 
	if (chopIdx == eigvals.size()){
		eigvals = eigvals.tail(1).eval();
		Z = Z.col(Z.cols()-1).eval();
	} else {
		int nLeftOver = eigvals.size()-chopIdx;
		eigvals = eigvals.tail(nLeftOver).eval();
		Z = Z.topRightCorner(Z.rows(), nLeftOver).eval();
	}
	//normalize the rows of Z
	const int nZCols = Z.cols(); //number of clusters currently instantiated
	for (int j = 0; j < nA; j++){
		double rownorm = sqrt(Z.row(j).squaredNorm());
		if (rownorm == 0){ //if rownorm is a hard zero, just set the row to ones -- 
							//it's the only thing we can do, lambda was set too high
			Z.row(j) = MXd::Ones(1, nZCols);
			rownorm = sqrt(Z.row(j).squaredNorm());
		}
		Z.row(j) *= (1.0/rownorm);
	}
	//initialize X (constrained version of Z) 
	MXd X(nA, nZCols);
	//initialize V (rotation matrix on Z to make Z*V close to X)
	MXd V(nZCols, nZCols);

	//propose nRestarts V trials
	V.setZero();
	//initialize unitary V via ``most orthogonal rows'' method
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> uniint(0, nA-1);
	int rndRow = uniint(gen);
	V.col(0) = Z.row(rndRow).transpose();
	MXd c(nA, 1);
	c.setZero();
	for (int j = 1; j < nZCols; j++){
		c += (Z*V.col(j-1)).cwiseAbs();
		int unused, nxtRow;
		c.minCoeff(&nxtRow, &unused);
		V.col(j) = Z.row(nxtRow).transpose();
	}
	this->orthonormalize(V);

	//solve the alternating minimization for X, V
	double obj, prevobj;
	obj = prevobj = numeric_limits<double>::infinity();
	do{
		prevobj = obj;
		MXd ZV = Z*V;
		//Solve for X via nonmaximum suppression
		X.setZero();
		for (int i = 0; i < nA; i++){
			int unused, maxCol;
			ZV.row(i).maxCoeff(&unused, &maxCol);
			X(i, maxCol) = 1;
		}

		//solve for V via SVD orthonormalization
		V = X.transpose()*Z;
		this->orthonormalize(V);

		obj = (X-Z*V).squaredNorm();
	} while( fabs(obj - prevobj)/obj > 1e-6);

	//Done -- pick out labels from X
	std::vector<int> lbls;
	for (int i = 0; i < nA; i++){
		int unused, maxCol;
		X.row(i).maxCoeff(&unused, &maxCol);
		lbls.push_back(maxCol);
	}
	return lbls;
}

template <typename G>
void KernDynMeans<G>::orthonormalize(MXd& V) const{ 
	//do a safe svd, with guaranteed orthonormal columns even in the event of numerically small singular values
	Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::HouseholderQRPreconditioner> svd(V, Eigen::ComputeFullU | Eigen::ComputeFullV); //slowest/safest preconditioning
	MXd U = svd.matrixU();
	MXd W = svd.matrixV();
	MXd id = MXd::Identity(V.rows(), V.cols());
	if ( (U.transpose()*U - id).squaredNorm() > 1e-6 || 
		 (W.transpose()*W - id).squaredNorm() > 1e-6  ){
		//buggy eigen sometimes spits out non-unitary U/W matrices
		//however, this only happens when the singular value is numerical-precision small, like 1e-34ish
		//so just pull those columns cause they don't matter anyway

		//find good columns
		vector<int> goodCols;
		for (int i = 0; i < U.cols(); i++){
			if ( fabs(U.col(i).squaredNorm()-1.0) < 1e-6 && fabs(W.col(i).squaredNorm()-1.0) < 1e-6){
				goodCols.push_back(i);
			}
		}
		//shift them up -- aliasing doesn't occur since goodcols increases by at least 1 each time
		for (int i = 0; i < goodCols.size(); i++){
			U.col(i) = U.col(goodCols[i]);
			W.col(i) = W.col(goodCols[i]);
		}
		//gram schmidt them just to be sure
		for (int i = 0; i < goodCols.size(); i++){
			for (int j = 0; j < i; j++){
				double dd = U.col(i).transpose()*U.col(j);
				U.col(i) -= U.col(j)*dd;
			}
			U.col(i) /= U.col(i).norm();
			for (int j = 0; j < i; j++){
				double dd = W.col(i).transpose()*W.col(j);
				W.col(i) -= W.col(j)*dd;
			}
			W.col(i) /= W.col(i).norm();
		}
		
		//graham schmidt to create the last few columns
		//U first
		int nGoodCols = goodCols.size();
		for (int i = 0; i < U.cols(); i++){
			if (nGoodCols == U.cols()){
				break;
			}
			VXd gscol = VXd::Zero(U.rows());
			gscol(i) = 1.0;
			//subtract the projections
			for (int j = 0; j < nGoodCols; j++){
				double dd = gscol.transpose()*U.col(j);
				gscol -= dd*U.col(j);
			}
			if (gscol.squaredNorm() > 1e-6){
				U.col(nGoodCols) = gscol/gscol.norm();
				nGoodCols++;
			}
		}
		if(nGoodCols != U.cols()){
			cout << "libkerndynmeans: ERROR: GRAHAM SCHMIDT DID NOT ORTHONORMALIZE U!" << endl;
			int ddd;
			cin >> ddd;
		}
		//now W
		nGoodCols = goodCols.size();
		for (int i = 0; i < W.cols(); i++){
			if (nGoodCols == W.cols()){
				break;
			}
			VXd gscol = VXd::Zero(U.rows());
			gscol(i) = 1.0;
			//subtract the projections
			for (int j = 0; j < nGoodCols; j++){
				double dd = gscol.transpose()*W.col(j);
				gscol -= dd*W.col(j);
			}
			if (gscol.squaredNorm() > 1e-6){
				W.col(nGoodCols) = gscol/gscol.norm();
				nGoodCols++;
			}
		}
		if(nGoodCols != W.cols()){
			cout << "libkerndynmeans: ERROR: GRAHAM SCHMIDT DID NOT ORTHONORMALIZE W!" << endl;
			int ddd;
			cin >> ddd;
		}
	}
	V  = W*U.transpose();
	return; 
}

template <typename G>
void KernDynMeans<G>::initializeSigma(const G& aff){
	int nNodes = aff.getNNodes();
	this->sigmaLB = 0.0;
	this->sigma = 0.0;
	this->sigmaUB = 0.0;
	int nOldPrms = aff.getNOldPrms();
	//check data rows
	for (int i = 0; i < nNodes; i++){
		double d = fabs(aff.diagSelfSimDD(i));
		double odsum = 0.0;
		for (int j = 0; j < i; j++){
			odsum += fabs(aff.simDD(i, j));
		}
		for (int j = i+1; j < nNodes; j++){
			odsum += fabs(aff.simDD(i, j));
		}
		for (int j = 0; j < nOldPrms; j++){
			odsum += fabs(aff.simDP(i, j));
		}
		if (odsum - d > this->sigmaUB){
			this->sigmaUB = odsum-d;
		}
	}
	//check old parameter rows
	for (int i = 0; i < nOldPrms; i++){
		double d = fabs(aff.selfSimPP(i));
		double odsum = 0.0;
		for (int j = 0; j < nNodes; j++){
			odsum += fabs(aff.simDP(j, i));
		}
		if ( this->gammas[i]*(odsum - d) > this->sigmaUB){
			this->sigmaUB = this->gammas[i]*(odsum - d);
		}
	}
	return;
}

template <class G>
template <typename T> void CoarseGraph<G>::coarsify(const T& aff){
	this->nOldPrms = aff.getNOldPrms();
	//coarsify and create the refinementmap
	//Pick a random order to traverse the data
	int nNodes = aff.getNNodes();
	std::vector<int> idcs(nNodes);
	std::iota(idcs.begin(), idcs.end(), 0);
	std::random_shuffle(idcs.begin(), idcs.end());
	//set up the vector to save which vertices have been marked
	std::vector<bool> marks(nNodes, false);
	this->refineMap.clear();
	for (int i = 0; i < idcs.size(); i++){
		int idxi = idcs[i];
		if (!marks[idxi]){//if the vertex hasn't already been merged to another
			double maxSim = 0;
			int maxId = -1;
			for (int j = i+1; j < idcs.size(); j++){//search all vertices after i (since all beforehave been merged)
				int idxj = idcs[j];
				if (!marks[idxj]){//only check it if it hasn't been marked
					double sim = aff.simDD(idxi, idxj);
					if (sim > maxSim && fabs(sim) > 1e-16){//1e-16 for keeping sparsity
						maxSim = sim;
						maxId = idxj;
					}
				}
			}
			//if maxId is still -1, then pair(i, -1) states correctly that i is a singleton
			this->refineMap.push_back(std::pair<int, int>(idxi, maxId));
			marks[idxi] = true;
			if (maxId >= 0){
				marks[maxId] = true;
			}
		}
	}
	int nNodesNew = this->refineMap.size();
	//now all merges have been found
	//create the coarsified graph -- only store upper triangle
	this->affdd = SMXd(nNodesNew, nNodesNew);
	this->affdp = SMXd(nNodesNew, this->nOldPrms);
	this->daffdd = VXd::Zero(nNodesNew);
	this->odaffdd = VXd::Zero(nNodesNew);
	this->affpp = VXd::Zero(this->nOldPrms);
	this->nodeCts = std::vector<int>(nNodesNew, 0);

	std::vector<TD> ddtrips, dptrips;
	for (int i = 0; i < this->refineMap.size(); i++){
		const int& i1 = this->refineMap[i].first;
		const int& i2 = this->refineMap[i].second;

		//sum up the node counts
		this->nodeCts[i] = aff.getNodeCt(i1) + (i2 != -1 ? aff.getNodeCt(i2) : 0);

		//off diagonal self similarity
		this->odaffdd(i) = aff.offDiagSelfSimDD(i1) + (i2 != -1 ? aff.offDiagSelfSimDD(i2) + aff.simDD(i1, i2) : 0);

		//diagonal self similarity
		this->daffdd(i) = aff.diagSelfSimDD(i1) + (i2 != -1 ? aff.diagSelfSimDD(i2) : 0);

		//data->data similarities
		for (int j = i+1; j < this->refineMap.size(); j++){
			const int& j1 = this->refineMap[j].first;
			const int& j2 = this->refineMap[j].second;
			double sim = aff.simDD(i1, j1) + (i2 != -1 ? aff.simDD(i2, j1) : 0) + (j2 != -1 ? aff.simDD(i1, j2) + (i2 != -1 ? aff.simDD(i2, j2) : 0) : 0);
			if (fabs(sim) > 1e-16){
				ddtrips.push_back(TD(i, j, sim));
			}
		}
		//data->param similarities
		for (int j = 0; j < this->nOldPrms; j++){
			double sim = aff.simDP(i1, j) + (i2 != -1 ? aff.simDP(i2, j) : 0);
			if (fabs(sim) > 1e-16){dptrips.push_back(TD(i, j, sim));}
		}
	}
	//param->param similarities
	for (int i = 0; i < this->nOldPrms; i++){
		this->affpp(i) = aff.selfSimPP(i);
	}

	//set the sparse matrices from triplets
	this->affdd.setFromTriplets(ddtrips.begin(), ddtrips.end());
	if (this->nOldPrms > 0){
		this->affdp.setFromTriplets(dptrips.begin(), dptrips.end());
	}
	//cout << "AFFDD: " << endl;
	//cout << MXd(this->affdd) << endl;
	//cout << "DAFFDD: " << endl;
	//cout << this->daffdd.transpose() << endl;
	//cout << "ODAFFDD: " << endl;
	//cout << this->odaffdd.transpose() << endl;
	//if (this->nOldPrms > 0){
	//	cout << "AFFDP: " << endl;
	//	cout << MXd(this->affdp) << endl;
	//	cout << "AFFPP: " << endl;
	//	cout << this->affpp.transpose() << endl;
	//}
	return;
}

template <class G>
double CoarseGraph<G>::diagSelfSimDD(const int i) const{
	return this->daffdd(i);
}
template <class G>
double CoarseGraph<G>::offDiagSelfSimDD(const int i) const{
	return this->odaffdd(i);
}
template <class G>
double CoarseGraph<G>::simDD(const int i, const int j) const{
	if(i == j){
		cout << "libkerndynmeans: ERROR: Do not use CoarseGraph::sim on indices i==j" << endl;
		cout << "libkerndynmeans: ERROR: Need to specify whether linear/quadratic self similarity." << endl;
		return 0.0;
	}
	if (i > j){//only upper triangle is stored
		return this->affdd.coeff(j, i);
	} else {
		return this->affdd.coeff(i, j);
	}
}

template <class G>
double CoarseGraph<G>::selfSimPP(const int i) const{
	return this->affpp(i);
}

template <class G>
double CoarseGraph<G>::simDP(const int i, const int j) const{
	return this->affdp.coeff(i, j);
}

template <class G>
std::vector<int> CoarseGraph<G>::getRefinedLabels(const std::vector<int>& lbls) const{
	//find the max index in merges to see how big the new labels should be
	int idxmax = -1;
	for (int i = 0; i < this->refineMap.size(); i++){
		if (this->refineMap[i].first > idxmax){
			idxmax = this->refineMap[i].first;
		}
		if (this->refineMap[i].second > idxmax){
			idxmax = this->refineMap[i].second;
		}
	}
	//fill in the extended labels by assigning all subnodes the label of the supernode
	std::vector<int> newlbls(idxmax+1, 0);
	for (int i = 0; i < lbls.size(); i++){
		newlbls[this->refineMap[i].first] = lbls[i];
		if (this->refineMap[i].second != -1){ //the -1 signal says that the node was just moved up a level in the hierarchy, no merge occured
										//so if there is a -1, just do nothing with merges.second
			newlbls[this->refineMap[i].second] = lbls[i];
		}
	}
	return newlbls;
}

template <class G>
int CoarseGraph<G>::getNodeCt(const int i) const{
	return this->nodeCts[i];
}

template <class G>
int CoarseGraph<G>::getNNodes() const{
	return this->affdd.cols();
}

template <class G>
int CoarseGraph<G>::getNOldPrms() const{
	return this->nOldPrms;
}

//std::vector<int> dynmLabelUpdate(std::vector<V2d> data, std::vector<int> lbls, std::vector<std::vector<int> > nodesets, double lambda){
//	vector<int> unqlbls = lbls;
//	sort(unqlbls.begin(), unqlbls.end());
//	unqlbls.erase(unique(unqlbls.begin(), unqlbls.end()), unqlbls.end());
//	int nextlbl = *max_element(lbls.begin(), lbls.end()) + 1;
//
//	std::map<int, V2d> prms;
//	for (int i = 0; i < unqlbls.size(); i++){
//		prms[unqlbls[i]] = V2d::Zero();
//		int cnt = 0;
//		for (int j = 0; j < lbls.size(); j++){
//			if (lbls[j] == unqlbls[i]){
//				for (int k = 0; k < nodesets[j].size(); k++){
//					prms[unqlbls[i]] += data[nodesets[j][k]];
//					cnt++;
//				}
//			}
//		}
//		prms[unqlbls[i]] /= cnt;
//	}
//	std::vector<int> newlbls(lbls.size(), 0);
//	for (int i = 0; i < lbls.size(); i++){
//		V2d newprm = V2d::Zero();
//		int cnt = 0;
//		for (int j = 0; j < nodesets[i].size(); j++){
//			newprm += data[nodesets[i][j]];
//			cnt++;
//		}
//		newprm /= cnt;
//		double minCost = lambda;
//		for (int j =0 ; j< nodesets[i].size(); j++){
//			minCost += (data[nodesets[i][j]]-newprm).squaredNorm();
//		}
//		int minLbl = -1;
//		for (int k = 0; k < unqlbls.size(); k++){
//			double cost = 0;
//			for (int j = 0; j < nodesets[i].size(); j++){
//				cost += (data[nodesets[i][j]]-prms[unqlbls[k]]).squaredNorm();
//			}
//			if (cost < minCost){
//				minCost = cost;
//				minLbl = unqlbls[k];
//			}
//		}
//		if (minLbl == -1){
//			minLbl = nextlbl;
//			nextlbl++;
//		}
//		//relabel the datapoint
//		newlbls[i] = minLbl;
//		//if a new cluster was created or old one revived, add stuff into unq/idcsinclus/inclustersum/nInClus/etc
//		if (find(unqlbls.begin(), unqlbls.end(), minLbl) == unqlbls.end()){
//			unqlbls.push_back(minLbl);
//			prms[minLbl] = newprm;
//		}
//	}
//	return newlbls;
//}
//
//void testLabelUpdate(){
//	double lambda = 0.05;
//	double T_Q = 6.8;
//	double K_tau = 1.01;
//	double Q = lambda/T_Q;
//	double tau = (T_Q*(K_tau-1.0)+1.0)/(T_Q-1.0);
//	int nRestarts = 10;
//	int nCoarsest = 10;
//	KernDynMeans<VectorGraph> kdm(lambda, Q, tau, true);
//
//	std::vector<V2d> data;
//	srand((unsigned)time(0));
//	std::default_random_engine gen;
//	std::normal_distribution<double> dist1(0.0, 0.03);
//	std::normal_distribution<double> dist2(1.0, 0.03);
//	for (int i = 0; i < 100; i++){
//		data.push_back(V2d::Zero());
//		data.back()(0) = dist1(gen);
//		data.back()(1) = dist1(gen);
//		data.push_back(V2d::Zero());
//		data.back()(0) = dist2(gen);
//		data.back()(1) = dist2(gen);
//	}
//	
//	VectorGraph vg;
//	vg.updateData(data);
//	CoarseGraph<VectorGraph> cg1, cg2, cg3;
//	cg1.coarsify(vg);
//	cg2.coarsify(cg1);
//	cg3.coarsify(cg2);
//
//	std::vector<int> lbls3 = kdm.baseCluster(cg3);
//	std::vector<int> lbls2 = cg3.getRefinedLabels(lbls3);
//	std::vector<int> lbls1 = cg2.getRefinedLabels(lbls2);
//	std::vector<int> lbls =  cg1.getRefinedLabels(lbls1);
//	std::vector<int> dmlbls3 = lbls3;
//	for (int i = 0; i < 10; i++){
//		cout << "Graph -- Labels: ";
//		for (int i = 0; i < lbls3.size(); i++){
//			cout << lbls3[i] << " ";
//		}
//		cout << "Obj: " << kdm.objective(cg3, lbls3) << endl;
//		lbls3 = kdm.updateLabels(cg3, lbls3);
//	}
//
//	std::vector< std::vector<int> > nodesets;
//	for (int i = 0; i < cg3.refineMap.size(); i++){
//		vector<int> idcs2, idcs1;
//		idcs2.push_back(cg3.refineMap[i].first);
//		if(cg3.refineMap[i].second != -1){idcs2.push_back(cg3.refineMap[i].second);}
//		for (int j =0; j < idcs2.size(); j++){
//			idcs1.push_back(cg2.refineMap[idcs2[j]].first);
//			if(cg2.refineMap[idcs2[j]].second != -1){idcs1.push_back(cg2.refineMap[idcs2[j]].second);}
//		}
//		nodesets.push_back(std::vector<int>());
//		for (int j =0; j < idcs1.size(); j++){
//			nodesets.back().push_back(cg1.refineMap[idcs1[j]].first);
//			if(cg1.refineMap[idcs1[j]].second != -1){nodesets.back().push_back(cg1.refineMap[idcs1[j]].second);}
//		}
//	}
//
//	for (int i = 0; i < 10; i++){
//		cout << "DynM -- Labels: ";
//		for (int i = 0; i < dmlbls3.size(); i++){
//			cout << dmlbls3[i] << " ";
//		}
//		cout << "Obj: " << kdm.objective(cg3, dmlbls3) << endl;
//		dmlbls3= dynmLabelUpdate(data, dmlbls3, nodesets, kdm.lambda);
//	}
//
//	exit(0);
//
//
//}

//template <typename G>
//template <typename T> std::vector<int> KernDynMeans<G>::clusterSplit(std::vector<T>& data, std::vector<int> lbls){
//	//pick a random cluster
//	vector<int> unqlbls = lbls;
//	sort(unqlbls.begin(), unqlbls.end());
//	unqlbls.erase(unique(unqlbls.begin(), unqlbls.end()), unqlbls.end());
//	std::random_device rd;
//	std::mt19937 gen(rd());
//	std::uniform_int_distribution<> randclus(0, unqlbls.size()-1);
//	int ksp = unqlbls[randclus(gen)];
//
//	//pick a random node within the cluster
//	std::vector<int> clusidcs;
//	for (int i = 0; i < lbls.size(); i++){
//		if (lbls[i] == ksp){
//			clusidcs.push_back(i);
//		}
//	}
//	if (clusidcs.size() == 1){
//		//no split, the cluster had only one node in it
//		return lbls;
//	}
//
//	double preobj = this->objective(data, lbls);
//
//	std::uniform_int_distribution<> randnode(0, clusidcs.size()-1);
//	int isp1 =randnode(gen);
//	int isp2 = isp1;
//	while(isp2 == isp1){
//		isp2 = randnode(gen);
//	}
//
//	//grow the two clusters by weighted greedy search
//	auto paircomp = []( std::pair<int, double> a, std::pair<int, double> b){return a.second > b.second;};
//	std::vector<int> c1, c2;
//	std::vector<double> maxSimTo(clusidcs.size(), -std::numeric_limits<double>::max());
//	std::vector<bool> maxSimToVia1(clusidcs.size(), true);
//	maxSimTo[isp1] = std::numeric_limits<double>::max();
//	maxSimtoVia1[isp1] = true;
//	maxSimto[isp2] = std::numeric_limits<double>::max();
//	maxSimtoVia1[isp2] = false;
//	while(c1.size()+c2.size() < clusidcs.size()){
//		//get the next best node to add to a cluster
//		int idx = std::distance(maxSimTo.begin(), std::max_element(maxSimTo.begin(), maxSimTo.end()));
//		//if it was linked from 1, add it to 1
//		if (maxSimToVia1[idx] == true){
//			c1.push_back(idx);
//			//update all the maxsims
//			for (int i = 0; i < maxSimTo.size(); i++){
//				double sim = data[clusidcs[idx]].sim(data[clusidcs[i]]);
//				if(sim > maxSimTo[i]){
//					maxSimTo[i] = sim;
//					maxSimToVia1[i] = true;
//				}
//			}
//			//set all the max sims to anything already in c1 to -inf, therefore they'll never be repicked
//			for (int i = 0; i < c1.size(); i++){
//				maxSimTo[c1[i]] = -std::numeric_limits<double>::max();
//			}
//		} else {
//			c2.push_back(idx);
//			//update all the maxsims
//			for (int i = 0; i < maxSimTo.size(); i++){
//				double sim = data[clusidcs[idx]].sim(data[clusidcs[i]]);
//				if(sim > maxSimTo[i]){
//					maxSimTo[i] = sim;
//					maxSimToVia1[i] = false;
//				}
//			}
//			//set all the max sims to anything already in c1 to -inf, therefore they'll never be repicked
//			for (int i = 0; i < c2.size(); i++){
//				maxSimTo[c2[i]] = -std::numeric_limits<double>::max();
//			}
//		}
//	}
//
//	//if the cluster was new, just create a new label for isp2's cluster 
//	if (std::find(this->oldprmlbls.begin(), this->oldprmlbls.end(), lbls[clusidcs[isp2]]) == this->oldprmlbls.end()){
//		std::vector<int> newlbls = lbls;
//		int oldlbl = lbls[clusidcs[isp2]];
//		int newlbl = this->nextlbl; 
//		for (int i = 0; i < c2.size(); i++){
//			newlbls[clusidcs[c2[i]]] = newlbl;
//		}
//		double postobj = this->objective(data, newlbls);
//		if (postobj < preobj){
//			this->nextlbl++;
//			return newlbls;
//		} else{
//			return lbls;
//		}
//	} else {
//		//if the cluster was linked to an old one, pick the least costly linkage
//		std::vector<int> newlbls1 = lbls;
//		std::vector<int> newlbls2 = lbls;
//		int oldlbl = lbls[clusidcs[isp2]];
//		int newlbl = this->nextlbl; 
//		for (int i = 0; i < c1.size(); i++){
//			newlbls2[clusidcs[c1[i]]] = newlbl;
//		}
//		for (int i = 0; i < c2.size(); i++){
//			newlbls2[clusidcs[c2[i]]] = newlbl;
//		}
//		double postobj1 = this->objective(data, newlbls1);
//		double postobj2 = this->objective(data, newlbls2);
//		if (postobj1 < postobj2 && postobj1 < preobj){
//			this->nextlbl++;
//			return newlbls1;
//		} else if (postobj2 < postobj1 && postobj2 < preobj){
//			this->nextlbl++;
//			return newlbls2;
//		} else {
//			return lbls;
//		}
//	}
//}
//
//template <typename G>
//template <typename T> std::vector<int> KernDynMeans<G>::clusterMerge(std::vector<T>& data, std::vector<int> lbls){
//
//	//pick two random clusters
//	vector<int> unqlbls = lbls;
//	sort(unqlbls.begin(), unqlbls.end());
//	unqlbls.erase(unique(unqlbls.begin(), unqlbls.end()), unqlbls.end());
//	if (unqlbls.size() == 1){
//		//no merge, only one cluster
//		return lbls;
//	}
//
//	double preobj = this->objective(data, lbls);
//
//	std::random_device rd;
//	std::mt19937 gen(rd());
//	std::uniform_int_distribution<> randclus(0, unqlbls.size()-1);
//	int km1 = unqlbls[randclus(gen)];
//	int km2 = km1;
//	while(km2 == km1){
//		km2 = unqlbls[randclus(gen)];
//	}
//	//try merging both ways (due to asymmetry when one/both clusters are related to old clusters)
//	std::vector<int> newlbls1 = lbls;
//	std::vector<int> newlbls2 = lbls;
//	for (int i = 0; i < newlbls1.size(); i++){
//		if (newlbls1[i] == km2){
//			newlbls1[i] = km1;
//		}
//	}
//	double postobj1 = this->objective(data, newlbls1);
//	for (int i = 0; i < newlbls2.size(); i++){
//		if (newlbls2[i] == km1){
//			newlbls2[i] = km2;
//		}
//	}
//	double postobj2 = this->objective(data, newlbls2);
//	if (postobj1 < postobj2 && postobj1 < preobj){
//		return newlbls1;
//	} else if (postobj2 < postobj1 && postobj2 < preobj){
//		return newlbls2;
//	} else{
//		return lbls;
//	}
//}


#define __KERNDYNMEANS_IMPL_HPP
#endif /* __KERNDYNMEANS_IMPL_HPP */
