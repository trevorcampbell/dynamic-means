#ifndef __KERNDYNMEANS_IMPL_HPP
template<typename D, typename C, typename P>
KernDynMeans<D,C,P>::KernDynMeans(double lambda, double Q, double tau, bool verbose){
	this->verbose = verbose;
	this->nextlbl = 0;
	this->ages.clear();
	this->oldprms.clear();
	this->oldprmlbls.clear();
	this->weights.clear();
	this->gammas.clear();
	this->agecosts.clear();
	this->lambda = lambda;
	this->Q = Q;
	this->tau = tau;
}

template<typename D, typename C, typename P>
KernDynMeans<D,C,P>::~KernDynMeans(){
}

template<typename D, typename C, typename P>
void KernDynMeans<D,C,P>::reset(){
	this->nextlbl = 0;
	this->ages.clear();
	this->oldprms.clear();
	this->oldprmlbls.clear();
	this->weights.clear();
	this->gammas.clear();
	this->agecosts.clear();
}

//This function updates the weights/ages of all the clusters after each clustering step is complete
template <typename D, typename P>
void KernDynMeans<D,C,P>::updateState(const vector<D>& data, const vector<int>& lbls){
	//first increment the age of everything
	//this will be undone below for any current clusters
	for (int i = 0; i < this->ages.size(); i++){
		this->ages[i]++;
	}
	//get a map from label -> vector of data
	map<int, vector<D> > m;
	for (int i = 0; i < data.size(); i++){
		m[lbls[i]].push_back(data[i]);
	}
	//for every current cluster
	for (auto it = m.begin(); it != m.end(); ++it){
		//check if it's an old label
		auto it2 = find(this->oldprmlbls.begin(), this->oldprmlbls.end(), it->first);
		if (it2 == this->oldprmlbls.end()){
			//it's a new cluster, so create stuff for it
			this->ages.push_back(1);
			this->weights.push_back(it->second.size());
			this->oldprmlbls.push_back(it->first);
			this->oldprms.push_back(P(it->second));
		} else {
			//it's an old cluster, so update old stuff
			int oldidx = distance(this->oldprmlbls.begin(), it2);
			this->ages[oldidx] = 1;
			this->weights[oldidx] = this->gammas[oldidx] + it->second.size();
			this->oldprms[oldidx].update(it->second, this->gammas[oldidx]);
		}
	}
	//update the Gammas
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
		if (this->agecosts[i] > this->lamb){
			this->weights.erase(this->weights.begin()+i);
			this->ages.erase(this->ages.begin()+i);
			this->gammas.erase(this->gammas.begin()+i);
			this->agecosts.erase(this->agecosts.begin()+i);
			this->oldprms.erase(this->oldprms.begin()+i);
			this->oldprmlbls.erase(this->oldprmlbls.begin()+i);
			i--;
		}
	}

	//update this->nextlbl if new clusters were created
	int maxlbl = *max_element(lbls.begin(), lbls.end());
	if (maxlbl >= this->nextlbl){
		this->nextlbl = maxlbl+1;
	}

	//done
	return;
}

template<typename D, typename C, typename P>
void KernDynMeans<D,C,P>::cluster(std::vector<D>& data, const int nRestarts, const int nCoarsest, std::vector<int>& finalLabels, double& finalObj, double& tTaken){
	timeval tStart;
	gettimeofday(&tStart, NULL);

	const int nB = this->oldprms.size();
	const int nA = data.size();

	if (data.size() <= 0){
		cout << "libkerndynmeans: WARNING: data size <=0 (= " << nA << "); Returning empty labels."<<  endl;
		finalLabels = vector<int>();
		timeTaken = 0;
		finalObj = 0;
		return;
	}
	if (nRestarts <= 0){
		cout << "libkerndynmeans: ERROR: nRestarts <=0 (= " << nRestarts << ")"<<  endl;
		return;
	}
	if (verbose){
		cout << "libkerndynmeans: Clustering " << nA << " datapoints with " << nRestarts << " restarts." << endl;
		cout << "libkerndynmeans: " << nB << " old clusters from previous timesteps." << endl;
	}
	std::vector<int> minLbls;
	double minObj = std::numeric_limits<double>::max();
	for(int rest = 0; rest < nRestarts; rest++){
		if (verbose){
			cout << "libkerndynmeans: Attempt " << rest+1 << "/" << nRestarts << ", Obj = " << minObj << endl;
			cout << "libkerndynmeans: Coarsifying " << nA << " nodes..." << endl;
		}

		//first, form the coarsification levels in the graph
		std::stack<std::vector<C> > coarsestack; //stores coarsified nodes
		std::stack<std::vector<std::pair<int, int> > > mergestack; //mergestack.top() stores the pairs that were merged to form coarsestack.top()
		auto crs = this->coarsify(data);
		coarsestack.push(crs.first);
		mergestack.push(crs.second);
		while(coarsestack.top().size() > nCoarsest){
			if (verbose){
				cout << "libkerndynmeans: Coarsifying " << coarsestack.top().size() << " nodes at level " << coarsestack.size() << "." << endl;
			}
			crs = this->coarsify(coarsestack.top());
			coarsestack.push(crs.first);
			mergestack.push(crs.second);
		}
		if (verbose){
			cout << "libkerndynmeans: Done coarsifying, top level " << coarsestack.size() << " has " << coarsestack.top().size() << " nodes." << endl;
		}
		//next, step down through the refinements and cluster, initializing from the coarser level
		std::vector<int> lbls;
		while(!coarsestack.empty()){
			if (verbose){
				cout << "libkerndynmeans: Running clustering at level " << coarsestack.size() << " with " << coarsestack.top.size() << " nodes." << endl;
			}
			//optimize the labels for the current top of coarsestack
			lbls = this->clusterAtLevel(coarsestack.top(), lbls); //lbls starts out empty, clusterAtLevel knows to use a base clustering
			coarsestack.pop();
			//distribute the labels to the next level down
			lbls = this->refine(mergestack.top(), lbls);
			mergestack.pop();
		}
		if (verbose){
			cout << "libkerndynmeans: Running final clustering at data level." << endl;
		}
		//final clustering at the data level
		lbls = this->clusterAtLevel(data, lbls);

		//finally, compute the kernelized dynamic means objective
		if (verbose){
			cout << "libkerndynmeans: Computing objective..." << endl;
		}
		double obj = this->objective(data, lbls);
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
	this->updateState(data, minLbls);

	//output results
	finalObj =  minObj;
	finalLabels = minLbls;
	//get final time taken
	timeval tCur;
	gettimeofday(&tCur, NULL);
	timeTaken = (double)(tCur.tv_sec - tStart.tv_sec) + (double)(tCur.tv_usec - tStart.tv_usec)/1.0e6;
	return;
}


template<typename D, typename C, typename P>
template <typename T> std::vector<int> KernDynMeans<D,C,P>::clusterAtLevel(std::vector<T>& data, std::vector<int> lbls){
	if (initlbls.size() < data.size()){ // Base Clustering -- Use spectral clustering on data, maximum bipartite matching to link old clusters
		//get the data labels from spectral clustering
		SpecDynMeans<T, P> sdm(this->lambda, this->Q, this->tau);
		double tmpobj = 0;
		double tmpt = 0;
 		sdm.cluster(data, 1, data.size(), SpecDynMeans::EigenSolverType::REDSVD, lbls, tmpobj, tmpt);

		//find the optimal correspondence between old/current clusters
		lbls = this->updateOldNewCorrespondence(data, lbls);
		//initlbls is now ready for regular refinement iterations
	}

	//run the refinement iterations
	double obj = this->objective(data, lbls);
	for (int i = 0; i < 500; i++){
		lbls = this->updateLabels(data, lbls);
		lbls = this->clusterSplit(data, lbls);
		lbls = this->clusterMerge(data, lbls);
		lbls = this->updateOldNewCorrespondence(data, lbls);
		obj = this->objective(data, lbls);
		if (verbose){ cout << "libkerndynmeans: Objective = " << obj << "                           \r" << flush;}
	}
	if (verbose){cout << endl;}
	return lbls;
}

template <typename D, typename C, typename P>
template <typename T> std::vector<int> KernDynMeans<D,C,P>::updateLabels(std::vector<T>& data, std::vector<int> lbls){
	//get the unique labels
	vector<int> unqlbls = lbls;
	sort(unqlbls.begin(), unqlbls.end());
	unqlbls.erase(unique(unqlbls.begin(), unqlbls.end()), unqlbls.end());

	//get the number of observations in each cluster
	std::map<int, std::vector<T> > dInClus;
	for (int i = 0; i < lbls.size(); i++){
		dInClus[lbls[i]].push_back(data[i]);
	}

	//minimize the cost associated with each observation individually based on the old labelling
	std::vector<int> newlbls(lbls.size(), 0);
	for (int i = 0; i < lbls.size(); i++){
		double minCost = std::numeric_limits<double>::max();
		int minLbl = -1;
		for (int k = 0; k < unqlbls.size(); k++){
			auto it = find(this->oldprmlbls.begin(), this->oldprmlbls.end(), unqlbls[k]);
			if (it == this->oldprmlbls.end()){
				//new cluster
				double cost = (1.0 - 1.0/dInClus[unqlbls[k]].size())*data[i].sim(data[i]);
				for (int j = 0; j < dInClus.size(); j++){
					cost += -2.0/dInClus[unqlbls[k]].size()*data[i].sim(data[j]);
				}
				if (cost < minCost){
					minCost = cost;
					minLbl = unqlbls[k];
				}
			} else {
				//old cluster
				int oldidx = std::distance(this->oldprmlbls.begin(), it);
				double cost = (1.0 - 1.0/(this->gammas[oldidx] + dInClus[unqlbls[k]].size()))*data[i].sim(data[i]);
				for (int j = 0; j < dInClus.size(); j++){
					cost += -2.0/(this->gammas[oldidx]+dInClus[unqlbls[k]].size())*data[i].sim(data[j]);
				}
				cost += -2.0*this->gammas[oldidx]/(this->gammas[oldidx]+dInClus[unqlbls[k]].size())*this->oldprms[oldidx].sim(data[i]);
				if (cost < minCost){
					minCost = cost;
					minLbl = unqlbls[k];
				}
			}
		}
		newlbls[i] = minLbl;
	}
	return newlbls;
}


template <typename D, typename C, typename P>
template <typename T> std::vector<int> KernDynMeans<D,C,P>::clusterSplit(std::vector<T>& data, std::vector<int> lbls){
	//pick a random cluster
	vector<int> unqlbls = lbls;
	sort(unqlbls.begin(), unqlbls.end());
	unqlbls.erase(unique(unqlbls.begin(), unqlbls.end()), unqlbls.end());
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> randclus(0, unqlbls.size()-1);
	int ksp = unqlbls[randclus(gen)];

	//pick a random node within the cluster
	std::vector<int> clusidcs;
	for (int i = 0; i < lbls.size(); i++){
		if (lbls[i] == ksp){
			clusidcs.push_back(i);
		}
	}
	if (clusidcs.size() == 1){
		//no split, the cluster had only one node in it
		return lbls;
	}

	double preobj = this->objective(data, lbls);

	std::uniform_int_distribution<> randnode(0, clusidcs.size()-1);
	int isp1 =randnode(gen);
	int isp2 = isp1;
	while(isp2 == isp1){
		isp2 = randnode(gen);
	}

	//grow the two clusters by weighted greedy search
	auto paircomp = []( std::pair<int, double> a, std::pair<int, double> b){return a.second > b.second;};
	std::vector<int> c1, c2;
	std::vector<double> maxSimTo(clusidcs.size(), -std::numeric_limits<double>::max());
	std::vector<bool> maxSimToVia1(clusidcs.size(), true);
	maxSimTo[isp1] = std::numeric_limits<double>::max();
	maxSimtoVia1[isp1] = true;
	maxSimto[isp2] = std::numeric_limits<double>::max();
	maxSimtoVia1[isp2] = false;
	while(c1.size()+c2.size() < clusidcs.size()){
		//get the next best node to add to a cluster
		int idx = std::distance(maxSimTo.begin(), std::max_element(maxSimTo.begin(), maxSimTo.end()));
		//if it was linked from 1, add it to 1
		if (maxSimToVia1[idx] == true){
			c1.push_back(idx);
			//update all the maxsims
			for (int i = 0; i < maxSimTo.size(); i++){
				double sim = data[clusidcs[idx]].sim(data[clusidcs[i]]);
				if(sim > maxSimTo[i]){
					maxSimTo[i] = sim;
					maxSimToVia1[i] = true;
				}
			}
			//set all the max sims to anything already in c1 to -inf, therefore they'll never be repicked
			for (int i = 0; i < c1.size(); i++){
				maxSimTo[c1[i]] = -std::numeric_limits<double>::max();
			}
		} else {
			c2.push_back(idx);
			//update all the maxsims
			for (int i = 0; i < maxSimTo.size(); i++){
				double sim = data[clusidcs[idx]].sim(data[clusidcs[i]]);
				if(sim > maxSimTo[i]){
					maxSimTo[i] = sim;
					maxSimToVia1[i] = false;
				}
			}
			//set all the max sims to anything already in c1 to -inf, therefore they'll never be repicked
			for (int i = 0; i < c2.size(); i++){
				maxSimTo[c2[i]] = -std::numeric_limits<double>::max();
			}
		}
	}

	//if the cluster was new, just create a new label for isp2's cluster 
	if (std::find(this->oldprmlbls.begin(), this->oldprmlbls.end(), lbls[clusidcs[isp2]]) == this->oldprmlbls.end()){
		std::vector<int> newlbls = lbls;
		int oldlbl = lbls[clusidcs[isp2]];
		int newlbl = this->nextlbl; 
		for (int i = 0; i < c2.size(); i++){
			newlbls[clusidcs[c2[i]]] = newlbl;
		}
		double postobj = this->objective(data, newlbls);
		if (postobj < preobj){
			this->nextlbl++;
			return newlbls;
		} else{
			return lbls;
		}
	} else {
		//if the cluster was linked to an old one, pick the least costly linkage
		std::vector<int> newlbls1 = lbls;
		std::vector<int> newlbls2 = lbls;
		int oldlbl = lbls[clusidcs[isp2]];
		int newlbl = this->nextlbl; 
		for (int i = 0; i < c1.size(); i++){
			newlbls2[clusidcs[c1[i]]] = newlbl;
		}
		for (int i = 0; i < c2.size(); i++){
			newlbls2[clusidcs[c2[i]]] = newlbl;
		}
		double postobj1 = this->objective(data, newlbls1);
		double postobj2 = this->objective(data, newlbls2);
		if (postobj1 < postobj2 && postobj1 < preobj){
			this->nextlbl++;
			return newlbls1;
		} else if (postobj2 < postobj1 && postobj2 < preobj){
			this->nextlbl++;
			return newlbls2;
		} else {
			return lbls;
		}
	}
}

template <typename D, typename C, typename P>
template <typename T> std::vector<int> KernDynMeans<D,C,P>::clusterMerge(std::vector<T>& data, std::vector<int> lbls){

	//pick two random clusters
	vector<int> unqlbls = lbls;
	sort(unqlbls.begin(), unqlbls.end());
	unqlbls.erase(unique(unqlbls.begin(), unqlbls.end()), unqlbls.end());
	if (unqlbls.size() == 1){
		//no merge, only one cluster
		return lbls;
	}

	double preobj = this->objective(data, lbls);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> randclus(0, unqlbls.size()-1);
	int km1 = unqlbls[randclus(gen)];
	int km2 = km1;
	while(km2 == km1){
		km2 = unqlbls[randclus(gen)];
	}
	//try merging both ways (due to asymmetry when one/both clusters are related to old clusters)
	std::vector<int> newlbls1 = lbls;
	std::vector<int> newlbls2 = lbls;
	for (int i = 0; i < newlbls1.size(); i++){
		if (newlbls1[i] == km2){
			newlbls1[i] = km1;
		}
	}
	double postobj1 = this->objective(data, newlbls1);
	for (int i = 0; i < newlbls2.size(); i++){
		if (newlbls2[i] == km1){
			newlbls2[i] = km2;
		}
	}
	double postobj2 = this->objective(data, newlbls2);
	if (postobj1 < postobj2 && postobj1 < preobj){
		return newlbls1;
	} else if (postobj2 < postobj1 && postobj2 < preobj){
		return newlbls2;
	} else{
		return lbls;
	}
}

template <typename D, typename C, typename P>
template <typename T> std::vector<int> KernDynMeans<D,C,P>::updateOldNewCorrespondence(std::vector<T>& data, std::vector<int> lbls){
	//get the unique labels
	vector<int> unqlbls = lbls;
	sort(unqlbls.begin(), unqlbls.end());
	unqlbls.erase(unique(unqlbls.begin(), unqlbls.end()), unqlbls.end());

	//get the number of observations in each cluster
	std::map<int, int> nInClus;
	for (int i = 0; i < unqlbls.size(); i++){
		nInClus[unqlbls[i]] = 0;
	}
	for (int i = 0; i < lbls.size(); i++){
		nInClus[lbls[i]]++;
	}


	//get the old/new correspondences from bipartite matching
 	vector< pair<int, int> > nodePairs; //new clusters in index 0, old clusters + one null cluster in index 1
 	vector< double > edgeWeights;
	for (int i = 0; i < unqlbls.size(); i++){
		for (int j = 0; j < this->oldprmlbls.size(); j++){
			nodePairs.push_back(std::pair<int, int>(unqlbls[i], this->oldprmlbls[j]) );
			double ewt = this->agecosts[j];
			ewt += this->gammas[j]*nInClus[unqlbls[i]]/(this->gammas[j]+nInClus[unqlbls[i]])*this->oldprms[j].sim(this->oldprms[j]);
			for (int k = 0; k < initlbls.size(); k++){
				if (lbls[k] == unqlbls[i]){
					ewt += -2.0*this->gammas[j]/(this->gammas[j]+nInClus[unqlbls[i]])*this->oldprms[j].sim(data[k]);
				}
			}
			edgeWeights.push_back(ewt);
		}
		//-1 is the new cluster option
		nodePairs.push_back( std::pair<int, int>(unqlbls[i], -1) );
		edgeWeights.push_back(this->lambda);
	}
	map<int, int> matching = this->getMinWtMatching(nodePairs, edgeWeights);

	//relabel initlbls based on the old/new correspondences
	for (auto it = matching.begin(); it != matching.end(); ++it){
		if (it->second != -1){ //if the current cluster isn't new
			//replace all labels in initlbls to the old cluster label
			for (int i = 0; i < initlbls.size(); i++){
				if (lbls[i] == it->first){
					lbls[i] = it->second;
				}
			}
		}
	}
	return lbls;
}

template <typename D, typename C, typename P>
map<int, int> KernDynMeans<D,C,P>::getMinWtMatching(vector< pair<int, int> > nodePairs, vector<double> edgeWeights ) const{
	//get params
	int nVars = edgeWeights.size();

	//start up GRB
	GRBEnv grbenv;
	grbenv.set(GRB_IntParam_OutputFlag, 0);//controls the output of Gurobi - 0 means no output, 1 means normal output
	grbenv.set(GRB_IntParam_Threads, 1);//controls the number of threads Gurobi uses - I force it to use 1 since the optimization in this algorithm is fairly small/simple
									     	//											and it just ends up wasting time constantly creating/deleting threads otherwise
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
	for (int i = 0; i < B.size(); i++){
		GRBLinExpr constrlhs;
		for (int j = 0; j < nVars; j++){
			if (nodePairs[j].second == B[i]){
				constrlhs += 1.0*grbvars[j];
			}
		}
		grbmodel.addConstr(constrlhs, GRB_LESS_EQUAL, 1);
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
		cout << e.getMessage() << endl;
	} catch (...){
		cout << "libkerndynmeans: ERROR: Unhandled Gurobi exception during optimization" << endl;
	}
}


template<typename D, typename C, typename P>
std::vector<int> KernDynMeans<D,C,P>::refine(std::vector< std::pair<int, int> > merges, std::vector<int> lbls){
	//find the max index in merges to see how big the new labels should be
	int lblmax = -1;
	for (auto it = merges.begin(); it != merges.end(); ++it){
		if (it->first > lblmax){
			lblmax = it->first;
		}
		if (it->second > lblmax){
			lblmax = it->second;
		}
	}
	//fill in the extended labels by assigning all subnodes the label of the supernode
	std::vector<int> newlbls(lblmax+1, 0);
	for (int i = 0; i < merges.size(); i++){
		newlbls[merges[i].first] = lbls[i];
		newlbls[merges[i].second] = lbls[i];
	}
	return newlbls;
}

template<typename D, typename C, typename P>
template<typename T> std::pair< std::vector<C>, std::vector<std::pair<int, int> > >  KernDynMeans<D,C,P>::coarsify(std::vector<T>& data){
	//Pick a random order to traverse the data
	std::vector<int> idcs(data.size());
	std::iota(idcs.begin(), idcs.end(), 0);
	std::random_shuffle(idcs.begin(), idcs.end());
	//set up the vector to save which vertices have been marked
	std::vector<bool> marks(data.size(), false);
	std::vector< std::pair<int, int> > merges;
	for (int i = 0; i < idcs.size(); i++){
		if (!marks[i]){//if the vertex hasn't already been merged to another
			double maxSim = 0;
			int maxId = -1;
			for (int j = i+1; j < idcs.size(); j++){//search all vertices after i (since all beforehave been merged)
				if (!marks[j]){
					double sim = data[i].sim(data[j]);
					if (sim > maxSim && sim > 1e-16){//1e-16 for keeping sparsity, only choose those whose mark is false
						maxSim = sim;
						maxId = j;
					}
				}
			}
			//if maxId is still -1, then pair(i, -1) states correctly that i is a singleton
			merges.push_back( std::pair<int, int>(i, maxId));
			marks[i] = true;
			if (maxId >= 0){
				marks[maxId] = true;
			}
		}
	}
	//now all merges have been created
	//create the coarsified nodes
	std::vector<C> coarse;
	for (int i = 0; i < merges.size(); i++){
		coarse.push_back( C(data[merges[i].first], data[merges[i].second]));
	}
	return std::pair< std::vector<C>, std::vector<std::pair<int, int> > >(coarse, merges);
}

template<typename D, typename C, typename P>
double KernDynMeans<D,C,P>::objective(std::vector<D>& data, std::vector<int> lbls){
	double cost = 0;
	//get a map from label to clusters
	map<int, vector<D> > m;
	for (int i = 0; i < data.size(); i++){
		m[lbls[i]].push_back(data[i]);
	}
	//for every current cluster
	for (auto it = m.begin(); it != m.end(); ++it){
		//check if it's an old label
		auto it2 = find(this->oldprmlbls.begin(), this->oldprmlbls.end(), it->first);
		if (it2 == this->oldprmlbls.end()){ //it's a new cluster
			cost += this->lambda;//new cluster penalty
			//ratio association term
			for (int i = 0; i < it->second.size(); i++){
				cost += (1.0 - 1.0/it->second.size())*it->second[i].sim(it->second[i]) //diagonal elements
				for (int j = i+1; j < it->second.size(); j++){
					cost -= 2.0/it->second.size()*it->second[i].sim(it->second[j]); //off-diagonal elements
				}
			}
		} else { //it's an old cluster
			int oldidx = distance(this->oldprmlbls.begin(), it2);
			cost += this->agecosts[oldidx];//old cluster penalty
			//ratio association term
			for (int i = 0; i < it->second.size(); i++){
				cost += (1.0 - 1.0/it->second.size())*it->second[i].sim(it->second[i]) //diagonal elements
				for (int j = i+1; j < it->second.size(); j++){
					cost += -2.0/it->second.size()*it->second[i].sim(it->second[j]); //off-diagonal elements
				}
			}
			cost += this->gammas[oldidx]*it->second.size()/(this->gammas[oldidx]+it->second.size())*this->oldprms[oldidx].sim(this->oldprms[oldidx]);//old prm self-similarity
			//old prm ratio association term
			for (int i = 0; i < it->second.size(); i++){
				cost += -2.0*this->gammas[oldidx]/(this->gammas[oldidx]+it->second.size())*this->oldprms[oldidx].sim(it->second[i]);
			}
		}
	}
	return cost;
}

#define __KERNDYNMEANS_IMPL_HPP
#endif /* __KERNDYNMEANS_IMPL_HPP */
