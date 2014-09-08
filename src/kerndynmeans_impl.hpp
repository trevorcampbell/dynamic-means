#ifndef __KERNDYNMEANS_IMPL_HPP


//TODO
//TODO
//TODO REMOVE BELOW THIS LINE UNTIL THE "TODO REMOVE ABOVE THIS LINE" WHEN DONE TESTING WITH DMEANS OBJ
//TODO
//TODO

typedef Eigen::Vector2d V2d;
template<typename D, typename C, typename P>
double KernDynMeans<D,C,P>::computedynmeansobj(std::vector<D>& data, std::vector<int> lbls){
		double obj = 0;
		vector<int> unqlbls = lbls;
		sort(unqlbls.begin(), unqlbls.end());
		unqlbls.erase(unique(unqlbls.begin(), unqlbls.end()), unqlbls.end());

		for (int i =0 ; i < unqlbls.size(); i++){
			V2d prm;
			V2d dmean;
			int n = 0;
			dmean.setZero();
			for (int j = 0; j < lbls.size(); j++){
				if (lbls[j] == unqlbls[i]){
					n++;
					dmean += data[j].v;
				}
			}
			dmean /= (double)n;
			auto it = find(this->oldprmlbls.begin(), this->oldprmlbls.end(), unqlbls[i]);
			if (it == this->oldprmlbls.end()){
				obj += this->lambda;
				prm = dmean;
			} else {
				int oldidx = std::distance(this->oldprmlbls.begin(), it);
				prm =  (this->oldprms[oldidx].v*this->gammas[oldidx] + n*dmean)/((double)(n+this->gammas[oldidx]));
				obj += this->agecosts[oldidx];
				obj += this->gammas[oldidx]*(this->oldprms[oldidx].v-prm).squaredNorm();
			}
			for (int j = 0; j < lbls.size(); j++){
				if (lbls[j] == unqlbls[i]){
					obj += (prm - data[j].v).squaredNorm();
				}
			}
		}
		return obj;
}
template<typename D, typename C, typename P>
double KernDynMeans<D,C,P>::computedynmeansobj(std::vector<C>& data, std::vector<int> lbls){
		double obj = 0;
		vector<int> unqlbls = lbls;
		sort(unqlbls.begin(), unqlbls.end());
		unqlbls.erase(unique(unqlbls.begin(), unqlbls.end()), unqlbls.end());

		for (int i =0 ; i < unqlbls.size(); i++){
			V2d prm;
			V2d dmean;
			int n = 0;
			dmean.setZero();
			for (int j = 0; j < lbls.size(); j++){
				if (lbls[j] == unqlbls[i]){
					n+= data[j].nv;
					for (int k = 0; k < data[j].vs.size(); k++){
						dmean += data[j].vs[k];
					}
				}
			}
			dmean /= (double)n;
			auto it = find(this->oldprmlbls.begin(), this->oldprmlbls.end(), unqlbls[i]);
			if (it == this->oldprmlbls.end()){
				obj += this->lambda;
				prm = dmean;
			} else {
				int oldidx = std::distance(this->oldprmlbls.begin(), it);
				prm =  (this->oldprms[oldidx].v*this->gammas[oldidx] + n*dmean)/((double)(n+this->gammas[oldidx]));
				obj += this->agecosts[oldidx];
				obj += this->gammas[oldidx]*(this->oldprms[oldidx].v-prm).squaredNorm();
			}
			for (int j = 0; j < lbls.size(); j++){
				if (lbls[j] == unqlbls[i]){
					for (int k = 0; k < data[j].vs.size(); k++){
						obj += (prm - data[j].vs[k]).squaredNorm();
					}
				}
			}
		}
		return obj;
}

template<typename D, typename C, typename P>
double KernDynMeans<D,C,P>::computedynmeansobj2(std::vector<D>& data, std::vector<int> lbls){
	vector<int> unqlbls = lbls;
	sort(unqlbls.begin(), unqlbls.end());
	unqlbls.erase(unique(unqlbls.begin(), unqlbls.end()), unqlbls.end());
	map<int, vector<D> > dInClus;
	map<int, double> nInClus;
	for (int i = 0; i < data.size(); i++){
		dInClus[lbls[i]].push_back(data[i]);
		if(nInClus.count(lbls[i]) == 0){
			nInClus[lbls[i]] = 0;
		}
		nInClus[lbls[i]] += data[i].getN();
	}

	double objective = 0;
	for (int i = 0; i < unqlbls.size(); i++){
			//add cost for new clusters - lambda
			// or add cost for old clusters - Q
			auto it = find(this->oldprmlbls.begin(), this->oldprmlbls.end(), unqlbls[i]);
			if (it != this->oldprmlbls.end()){
				int oldidx = std::distance(this->oldprmlbls.begin(), it);
				objective += this->Q*this->ages[oldidx];
			} else {
				objective += this->lambda;
			}
			V2d tmpvec;
			tmpvec.setZero();
			for (int j = 0; j < dInClus[unqlbls[i]].size(); j++){
				tmpvec = tmpvec + dInClus[unqlbls[i]][j].v;
			}
			tmpvec = tmpvec / nInClus[unqlbls[i]];
			V2d prm;
			prm.setZero();
			if (it != this->oldprmlbls.end()){
				int oldidx =  std::distance(this->oldprmlbls.begin(), it);
				prm = (this->oldprms[oldidx].v*this->gammas[oldidx] + tmpvec*nInClus[unqlbls[i]])/(this->gammas[oldidx] + nInClus[unqlbls[i]]);
				//add parameter lag cost
				double tmpsqdist = (prm - this->oldprms[oldidx].v).squaredNorm();
				objective += this->gammas[oldidx]*tmpsqdist;
			} else { //just setting a new param
				prm = tmpvec;
				//no lag cost for new params
			}
			//get cost for prms[i]
			for (int j = 0; j < dInClus[unqlbls[i]].size(); j++){
				objective += (prm - dInClus[unqlbls[i]][j].v).squaredNorm();
			}
	}
	return objective;
}

template<typename D, typename C, typename P>
double KernDynMeans<D,C,P>::computedynmeansobj2(std::vector<C>& data, std::vector<int> lbls){
	vector<int> unqlbls = lbls;
	sort(unqlbls.begin(), unqlbls.end());
	unqlbls.erase(unique(unqlbls.begin(), unqlbls.end()), unqlbls.end());
	map<int, vector<C> > dInClus;
	map<int, double> nInClus;
	for (int i = 0; i < data.size(); i++){
		dInClus[lbls[i]].push_back(data[i]);
		if(nInClus.count(lbls[i]) == 0){
			nInClus[lbls[i]] = 0;
		}
		nInClus[lbls[i]] += data[i].getN();
	}

	double objective = 0;
	for (int i = 0; i < unqlbls.size(); i++){
			//add cost for new clusters - lambda
			// or add cost for old clusters - Q
			if (find(this->oldprmlbls.begin(), this->oldprmlbls.end(), unqlbls[i]) != this->oldprmlbls.end()){
				objective += this->Q*this->ages[i];
			} else {
				objective += this->lambda;
			}
			V2d tmpvec;
			tmpvec.setZero();
			for (int j = 0; j < dInClus[unqlbls[i]].size(); j++){
				for (int k = 0; k < dInClus[unqlbls[i]][j].vs.size(); k++){
					tmpvec = tmpvec + dInClus[unqlbls[i]][j].vs[k];
				}
			}
			tmpvec = tmpvec / nInClus[unqlbls[i]];
			V2d prm;
			prm.setZero();
			auto it = find(this->oldprmlbls.begin(), this->oldprmlbls.end(), unqlbls[i]);
			if (it != this->oldprmlbls.end()){
				int oldidx =  std::distance(this->oldprmlbls.begin(), it);
				prm = (this->oldprms[oldidx].v*this->gammas[oldidx] + tmpvec*nInClus[unqlbls[i]])/(this->gammas[oldidx] + nInClus[unqlbls[i]]);
				//add parameter lag cost
				double tmpsqdist = (prm - this->oldprms[oldidx].v).squaredNorm();
				objective += this->gammas[oldidx]*tmpsqdist;
			} else { //just setting a new param
				prm = tmpvec;
				//no lag cost for new params
			}
			//get cost for prms[i]
			for (int j = 0; j < dInClus[unqlbls[i]].size(); j++){
				for (int k = 0; k < dInClus[unqlbls[i]][j].vs.size(); k++){
					objective += (prm - dInClus[unqlbls[i]][j].vs[k]).squaredNorm();
				}
			}
	}
	return objective;
}



//TODO
//TODO
//TODO REMOVE ABOVE THIS LINE UNTIL THE "TODO REMOVE BELOW THIS LINE" WHEN DONE TESTING WITH DMEANS OBJ
//TODO
//TODO

template<typename D, typename C, typename P>
KernDynMeans<D,C,P>::KernDynMeans(double lambda, double Q, double tau, bool verbose){
	if (lambda < 0 || Q < 0 || tau < 0){
		cout << "libkerndynmeans: ERROR: Parameters of Kernel Dynamic Means cannot be < 0." << endl;
		cout << "libkerndynmeans: Lambda: " << lambda << " Q: " << Q << " tau: " << tau << endl;
	}
	this->verbose = verbose;
	this->maxLblPrevUsed = -1;
	this->ages.clear();
	this->oldprms.clear();
	this->oldprmlbls.clear();
	this->weights.clear();
	this->gammas.clear();
	this->agecosts.clear();
	this->lambda = lambda;
	this->Q = Q;
	this->tau = tau;
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

template<typename D, typename C, typename P>
KernDynMeans<D,C,P>::~KernDynMeans(){
	delete this->grbenv;
}

template<typename D, typename C, typename P>
void KernDynMeans<D,C,P>::reset(){
	this->maxLblPrevUsed = 0;
	this->ages.clear();
	this->oldprms.clear();
	this->oldprmlbls.clear();
	this->weights.clear();
	this->gammas.clear();
	this->agecosts.clear();
}

//This function updates the weights/ages of all the clusters after each clustering step is complete
template <typename D, typename C, typename P>
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
		if (this->agecosts[i] > this->lambda){
			this->weights.erase(this->weights.begin()+i);
			this->ages.erase(this->ages.begin()+i);
			this->gammas.erase(this->gammas.begin()+i);
			this->agecosts.erase(this->agecosts.begin()+i);
			this->oldprms.erase(this->oldprms.begin()+i);
			this->oldprmlbls.erase(this->oldprmlbls.begin()+i);
			i--;
		}
	}

	//update this->maxLblPrevUsed if new clusters were created
	int maxlbl = *max_element(lbls.begin(), lbls.end());
	if (maxlbl >= this->maxLblPrevUsed){
		this->maxLblPrevUsed = maxlbl;
	}

	//done
	return;
}

template<typename D, typename C, typename P>
void KernDynMeans<D,C,P>::cluster(std::vector<D>& data, const int nRestarts, const int nCoarsest, std::vector<int>& finalLabels, double& finalObj, double& tTaken){
	timeval tStart;
	gettimeofday(&tStart, NULL);

	if (data.size() <= 0){
		cout << "libkerndynmeans: WARNING: data size <=0 (= " << data.size() << "); Returning empty labels."<<  endl;
		finalLabels = vector<int>();
		tTaken = 0;
		finalObj = 0;
		return;
	}
	if (nRestarts <= 0){
		cout << "libkerndynmeans: ERROR: nRestarts <=0 (= " << nRestarts << ")"<<  endl;
		return;
	}
	if (verbose){
		cout << "libkerndynmeans: Clustering " << data.size() << " datapoints with " << nRestarts << " restarts." << endl;
		cout << "libkerndynmeans: " << this->oldprms.size() << " old clusters from previous timesteps." << endl;
	}
	std::vector<int> minLbls;
	double minObj = std::numeric_limits<double>::max();
	for(int rest = 0; rest < nRestarts; rest++){
		if (verbose){
			cout << "libkerndynmeans: Attempt " << rest+1 << "/" << nRestarts << ", Minimum Obj = " << minObj << endl;
		}

		//first, form the coarsification levels in the graph if necessary
		std::vector<int> lbls;
		if(data.size() > nCoarsest){
			if (verbose){
				cout << "libkerndynmeans: Coarsifying " << data.size() << " nodes at level 0." << endl;
			}
			std::stack<std::vector<C> > coarsestack; //stores coarsified nodes
			std::stack<std::vector<std::pair<int, int> > > mergestack; //mergestack.top() stores the pairs that were merged to form coarsestack.top()
			auto crs = this->coarsify(data);
			coarsestack.push(crs.first);
			mergestack.push(crs.second);
			while(coarsestack.top().size() > nCoarsest){
				if (verbose){
					cout << "libkerndynmeans: Coarsifying " << coarsestack.top().size() << " nodes at level " << coarsestack.size() << "." << endl;
				}
				auto crs = this->coarsify(coarsestack.top());
				coarsestack.push(crs.first);
				mergestack.push(crs.second);
			}
			if (verbose){
				cout << "libkerndynmeans: Done coarsifying, top level " << coarsestack.size() << " has " << coarsestack.top().size() << " nodes." << endl;
			}
			//next, step down through the refinements and cluster, initializing from the coarser level
			while(!coarsestack.empty()){
				if (verbose){
					cout << "libkerndynmeans: Running clustering at level " << coarsestack.size() << " with " << coarsestack.top().size() << " nodes." << endl;
				}
				//optimize the labels for the current top of coarsestack
				lbls = this->clusterAtLevel(coarsestack.top(), lbls); //lbls starts out empty, clusterAtLevel knows to use a base clustering
				coarsestack.pop();
				//distribute the labels to the next level down
				lbls = this->refine(mergestack.top(), lbls);
				mergestack.pop();
			}
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
	tTaken = (double)(tCur.tv_sec - tStart.tv_sec) + (double)(tCur.tv_usec - tStart.tv_usec)/1.0e6;
	return;
}


template<typename D, typename C, typename P>
template <typename T> 
std::vector<int> KernDynMeans<D,C,P>::clusterAtLevel(std::vector<T>& data, std::vector<int> lbls){
	if (lbls.size() < data.size()){ // Base Clustering -- Use spectral clustering on data, maximum bipartite matching to link old clusters
		if (verbose){ cout << "Running base spectral clustering..." << endl;}
		//get the data labels from spectral clustering
		lbls = this->baseCluster(data);
		//find the optimal correspondence between old/current clusters
		lbls = this->updateOldNewCorrespondence(data, lbls);
		//initlbls is now ready for regular refinement iterations
		if (verbose){ cout << "Done base spectral clustering with objective: " << this->objective(data, lbls) << endl;}
	}

	//run the refinement iterations
	double prevobj = this->objective(data, lbls);
	double diff = 1.0;
	int itr = 0;
	while(diff > 1e-6){
		itr++;
		cout << "prelbl obj: " << this->objective(data, lbls) << endl;
		cout << "prelbl dmobj: " << this->computedynmeansobj(data, lbls) << endl;
		cout << "prelbl dmobj2: " << this->computedynmeansobj2(data, lbls) << endl;
		lbls = this->updateLabels(data, lbls);
		cout << "postlbl, preoldnew obj: " << this->objective(data, lbls) << endl;
		cout << "postlbl, preoldnew dmobj: " << this->computedynmeansobj(data, lbls) << endl;
		cout << "postlbl, preoldnew dmobj2: " << this->computedynmeansobj2(data, lbls) << endl;
		lbls = this->updateOldNewCorrespondence(data, lbls);
		cout << "postoldnew obj: " << this->objective(data, lbls) << endl;
		cout << "postoldnew dmobj: " << this->computedynmeansobj(data, lbls) << endl;
		cout << "postoldnew dmobj2: " << this->computedynmeansobj2(data, lbls) << endl;
		double obj = this->objective(data, lbls);
		diff = fabs((obj-prevobj)/obj);
		prevobj = obj;
		if (verbose){ cout << "libkerndynmeans: Kernelized clustering iteration " << itr << ", obj = " << obj << endl;}
	}
	if (verbose){cout << endl;}
	return lbls;
}

template <typename D, typename C, typename P>
template <typename T> 
std::vector<int> KernDynMeans<D,C,P>::updateLabels(std::vector<T>& data, std::vector<int> lbls){
	//get the unique labels
	vector<int> unqlbls = lbls;
	sort(unqlbls.begin(), unqlbls.end());
	unqlbls.erase(unique(unqlbls.begin(), unqlbls.end()), unqlbls.end());


	//get the observations in each cluster
	//and the sizes of each cluster
	map<int, vector<T> > dInClus;
	map<int, double> nInClus;
	for (int i = 0; i < data.size(); i++){
		dInClus[lbls[i]].push_back(data[i]);
		if(nInClus.count(lbls[i]) == 0){
			nInClus[lbls[i]] = 0;
		}
		nInClus[lbls[i]] += data[i].getN();
	}

	//precompute the squared cluster sums
	//and the old cluster sums
	std::map<int, double> sqClusterSum;
	std::map<int, double> oldClusterSum;
	for (int i = 0; i < unqlbls.size(); i++){
		double sqsum = 0;
		const std::vector<T>& clus = dInClus[unqlbls[i]];
		const int& lbl = unqlbls[i];
		for (int k = 0; k < clus.size(); k++){
			sqsum += clus[k].sim(clus[k]);
			for (int m = k+1; m < clus.size(); m++){
				sqsum += 2.0*clus[k].sim(clus[m]);
			}
		}
		sqClusterSum[lbl] = sqsum;
		auto it = find(this->oldprmlbls.begin(), this->oldprmlbls.end(), lbl);
		if (it != this->oldprmlbls.end()){
			int oldidx = std::distance(this->oldprmlbls.begin(), it);	
			double oldsum = 0;
			for (int k = 0; k < clus.size(); k++){
				oldsum += this->oldprms[oldidx].sim(clus[k]);
			}
			oldClusterSum[lbl] = oldsum;
		}
	}

	//minimize the cost associated with each observation individually based on the old labelling
	std::vector<int> newlbls(lbls.size(), 0);
	int nextlbl = this->maxLblPrevUsed+1;//for this round, handles labelling of new clusters
	for (int i = 0; i < lbls.size(); i++){
		double minCost = this->lambda; //default to creating a new cluster, and then try to beat it 
		int minLbl = -1;
		for (int k = 0; k < unqlbls.size(); k++){
			const std::vector<T>& clus = dInClus[unqlbls[k]];
			const int& lbl = unqlbls[k];
			auto it = find(this->oldprmlbls.begin(), this->oldprmlbls.end(), lbl);
			if (it == this->oldprmlbls.end()){
				//new instantiated cluster
				double factor = 1.0/nInClus[lbl];
				double cost = data[i].sim(data[i]) + factor*factor*sqClusterSum[lbl]; 
				for (int j = 0; j < clus.size(); j++){
					cost += -2.0*factor*data[i].sim(clus[j]);
				}
				if (cost < minCost){
					minCost = cost;
					minLbl = lbl;
				}
			} else {
				//old instantiated cluster
				int oldidx = std::distance(this->oldprmlbls.begin(), it);
				double factor = 1.0/(this->gammas[oldidx] + nInClus[lbl]);
				double cost = data[i].sim(data[i]) + factor*factor*sqClusterSum[lbl];
				cost += 2.0*this->gammas[oldidx]*factor*factor*oldClusterSum[lbl];
				cost += this->gammas[oldidx]*this->gammas[oldidx]*factor*factor*this->oldprms[oldidx].sim(this->oldprms[oldidx]);
				for (int j = 0; j < clus.size(); j++){
					cost += -2.0*factor*data[i].sim(clus[j]);
				}
				cost += -2.0*this->gammas[oldidx]*factor*this->oldprms[oldidx].sim(data[i]);
				if (cost < minCost){
					minCost = cost;
					minLbl = lbl;
				}
			}
		}
		//run through old uninstantiated clusters
		for (int k = 0; k < this->oldprmlbls.size(); k++){
			auto it = find(unqlbls.begin(), unqlbls.end(), this->oldprmlbls[k]);
			if (it == unqlbls.end()){
				double cost = this->agecosts[k] + 
						this->gammas[k]/(this->gammas[k]+1.0)*(data[i].sim(data[i])-2.0*this->oldprms[k].sim(data[i])+this->oldprms[k].sim(this->oldprms[k]));
				if (cost < minCost){
					minCost = cost;
					minLbl = this->oldprmlbls[k];
				}
			}
		}
		if (minLbl == -1){
			//create a new cluster
			newlbls[i] = nextlbl;
			unqlbls.push_back(nextlbl);
			dInClus[nextlbl].push_back(data[i]); //no problem with duplicating the datapoint here (it will still exist in dInClus in the old cluster)
													//this is because the distance comparisons to previous cluster centers shouldn't be affected by creating new centers
													//ergo, leave the old dInClus alone, but create a new entry so that other observations can switch tothis cluster
			sqClusterSum[nextlbl] = data[i].sim(data[i]);
			nextlbl++;
		} else {
			newlbls[i] = minLbl;
		}
	}
	return newlbls;
}


template <typename D, typename C, typename P>
template <typename T> 
std::vector<int> KernDynMeans<D,C,P>::updateOldNewCorrespondence(std::vector<T>& data, std::vector<int> lbls){
	//get the unique labels
	vector<int> unqlbls = lbls;
	sort(unqlbls.begin(), unqlbls.end());
	unqlbls.erase(unique(unqlbls.begin(), unqlbls.end()), unqlbls.end());

	//get the observations in each cluster
	//and the sizes of each cluster
	map<int, vector<T> > dInClus;
	map<int, double> nInClus;
	for (int i = 0; i < data.size(); i++){
		dInClus[lbls[i]].push_back(data[i]);
		if(nInClus.count(lbls[i]) == 0){
			nInClus[lbls[i]] = 0;
		}
		nInClus[lbls[i]] += data[i].getN();
	}
	//compute the squared cluster sums
	//no point computing old cluster sums, will just have to compute below once for all clusters anyway
	std::map<int, double> sqClusterSum;
	for (int i = 0; i < unqlbls.size(); i++){
		double sqsum = 0;
		const std::vector<T>& clus = dInClus[unqlbls[i]];
		const int& lbl = unqlbls[i];
		for (int k = 0; k < clus.size(); k++){
			sqsum += clus[k].sim(clus[k]);
			for (int m = k+1; m < clus.size(); m++){
				sqsum += 2.0*clus[k].sim(clus[m]);
			}
		}
		sqClusterSum[lbl] = sqsum;
	}

	//get the old/new correspondences from bipartite matching
 	vector< pair<int, int> > nodePairs; //new clusters in .first, old clusters + one null cluster in .second
 	vector< double > edgeWeights;
	for (int i = 0; i < unqlbls.size(); i++){
		const std::vector<T>& clus = dInClus[unqlbls[i]];
		const int& lbl = unqlbls[i];
		for (int j = 0; j < this->oldprmlbls.size(); j++){
			double ewt = this->agecosts[j]
						+ this->gammas[j]*nInClus[lbl]/(this->gammas[j]+nInClus[lbl])*this->oldprms[j].sim(this->oldprms[j])
						-1.0/(this->gammas[j]+nInClus[lbl])*sqClusterSum[lbl];
			double oldClusSum = 0;
			for (int k = 0; k < clus.size(); k++){
				oldClusSum += this->oldprms[j].sim(clus[k]);
			}
			ewt += -2.0*this->gammas[j]/(this->gammas[j]+nInClus[lbl])*oldClusSum;
			nodePairs.push_back(std::pair<int, int>(lbl, this->oldprmlbls[j]) );
			edgeWeights.push_back(ewt);
		}
		//-1 is the new cluster option
		nodePairs.push_back( std::pair<int, int>(lbl, -1) );
		edgeWeights.push_back(this->lambda-1.0/nInClus[lbl]*sqClusterSum[lbl]);
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

template <typename D, typename C, typename P>
map<int, int> KernDynMeans<D,C,P>::getMinWtMatching(vector< pair<int, int> > nodePairs, vector<double> edgeWeights ) const{
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
		if (merges[i].second != -1){ //the -1 signal says that the node was just moved up a level in the hierarchy, no merge occured
										//so if there is a -1, just do nothing with merges.second
			newlbls[merges[i].second] = lbls[i];
		}
	}
	return newlbls;
}

template<typename D, typename C, typename P>
template<typename T> 
std::pair< std::vector<C>, std::vector<std::pair<int, int> > >  KernDynMeans<D,C,P>::coarsify(std::vector<T>& data){
	//Pick a random order to traverse the data
	std::vector<int> idcs(data.size());
	std::iota(idcs.begin(), idcs.end(), 0);
	std::random_shuffle(idcs.begin(), idcs.end());
	//set up the vector to save which vertices have been marked
	std::vector<bool> marks(data.size(), false);
	std::vector< std::pair<int, int> > merges;
	for (int i = 0; i < idcs.size(); i++){
		int idxi = idcs[i];
		if (!marks[idxi]){//if the vertex hasn't already been merged to another
			double maxSim = 0;
			int maxId = -1;
			for (int j = i+1; j < idcs.size(); j++){//search all vertices after i (since all beforehave been merged)
				int idxj = idcs[j];
				if (!marks[idxj]){//only check it if it hasn't been marked
					double sim = data[idxi].sim(data[idxj]);
					if (sim > maxSim && sim > 1e-16){//1e-16 for keeping sparsity
						maxSim = sim;
						maxId = idxj;
					}
				}
			}
			//if maxId is still -1, then pair(i, -1) states correctly that i is a singleton
			merges.push_back(std::pair<int, int>(idxi, maxId));
			marks[idxi] = true;
			if (maxId >= 0){
				marks[maxId] = true;
			}
		}
	}
	//now all merges have been created
	//create the coarsified nodes
	std::vector<C> coarse;
	for (int i = 0; i < merges.size(); i++){
		if (merges[i].second != -1){ //the only one that can be -1 is the second entry
			coarse.push_back( C(data[merges[i].first], data[merges[i].second]));
		} else {
			coarse.push_back( C(data[merges[i].first]) );
		}
	}
	return std::pair< std::vector<C>, std::vector<std::pair<int, int> > >(coarse, merges);
}

template<typename D, typename C, typename P>
template<typename T> 
double KernDynMeans<D,C,P>::objective(std::vector<T>& data, std::vector<int> lbls){
	double cost = 0;
	//get a map from label to clusters
	map<int, vector<T> > dInClus;
	map<int, double> nInClus;
	for (int i = 0; i < data.size(); i++){
		dInClus[lbls[i]].push_back(data[i]);
		if(nInClus.count(lbls[i]) == 0){
			nInClus[lbls[i]] = 0;
		}
		nInClus[lbls[i]] += data[i].getN();
	}
	//for every current cluster
	for (auto it = dInClus.begin(); it != dInClus.end(); ++it){
		//check if it's an old label
		auto it2 = find(this->oldprmlbls.begin(), this->oldprmlbls.end(), it->first);
		const int& lbl = it->first;
		const std::vector<T>& clus = it->second;
		if (it2 == this->oldprmlbls.end()){ //it's a new cluster
			cost += this->lambda;//new cluster penalty
			//ratio association term
			for (int i = 0; i < clus.size(); i++){
				cost += (1.0 - 1.0/nInClus[lbl])*clus[i].sim(clus[i]); //diagonal elements
				for (int j = i+1; j < clus.size(); j++){
					cost -= 2.0/nInClus[lbl]*clus[i].sim(clus[j]); //off-diagonal elements
				}
			}
		} else { //it's an old cluster
			int oldidx = distance(this->oldprmlbls.begin(), it2);
			cost += this->agecosts[oldidx];//old cluster penalty
			//ratio association term
			for (int i = 0; i < clus.size(); i++){
				cost += (1.0 - 1.0/(this->gammas[oldidx]+nInClus[lbl]))*clus[i].sim(clus[i]); //diagonal elements
				for (int j = i+1; j < clus.size(); j++){
					cost += -2.0/(this->gammas[oldidx]+nInClus[lbl])*clus[i].sim(clus[j]); //off-diagonal elements
				}
			}
			cost += this->gammas[oldidx]*nInClus[lbl]/(this->gammas[oldidx]+nInClus[lbl])*this->oldprms[oldidx].sim(this->oldprms[oldidx]);//old prm self-similarity
			//old prm ratio association term
			for (int i = 0; i < clus.size(); i++){
				cost += -2.0*this->gammas[oldidx]/(this->gammas[oldidx]+nInClus[lbl])*this->oldprms[oldidx].sim(clus[i]);
			}
		}
	}
	return cost;
}

template<typename D, typename C, typename P>
template<typename T>
std::vector<int> KernDynMeans<D,C,P>::baseCluster(std::vector<T>& data){
	//compute the kernel matrix
	int nA = data.size();
	MXd K = MXd(nA, nA);
	K.setZero();
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

template <typename D, typename C, typename P>
void KernDynMeans<D,C,P>::orthonormalize(MXd& V) const{ 
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
			cout << "libspecdynmeans: ERROR: GRAHAM SCHMIDT DID NOT ORTHONORMALIZE U!" << endl;
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
			cout << "libspecdynmeans: ERROR: GRAHAM SCHMIDT DID NOT ORTHONORMALIZE W!" << endl;
			int ddd;
			cin >> ddd;
		}
	}
	V  = W*U.transpose();
	return; 
}



template <class G>
CoarseGraph<G>::CoarseGraph(const G& aff){
	this->coarsify(aff);
	this->oldPrmLbls = aff.getOldPrmLbls();
}

template <class G>
CoarseGraph<G>::CoarseGraph(const CoarseGraph& aff){
	this->coarsify(aff);
	this->oldPrmLbls = aff.getOldPrmLbls();
}


template <class G>
template <typename T> void CoarseGraph<G>::coarsify(const T& aff){
	//coarsify and create the refinementmap
	//Pick a random order to traverse the data
	int nData = aff.nData();
	int nOldPrm = aff.nOldPrm();
	std::vector<int> idcs(nData);
	std::iota(idcs.begin(), idcs.end(), 0);
	std::random_shuffle(idcs.begin(), idcs.end());
	//set up the vector to save which vertices have been marked
	std::vector<bool> marks(nData, false);
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
					if (sim > maxSim && sim > 1e-16){//1e-16 for keeping sparsity
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
	int nDataNew = this->refineMap.size();
	//now all merges have been found
	//create the coarsified graph -- only store upper triangle
	this->affdd = SMXd(nDataNew, nDataNew);
	this->affdp = SMXd(nDataNew, nOldPrm);
	this->linaffdd = VXd::Zero(nDataNew);
	this->affpp = VXd::Zero(nOldPrm);

	std::vector<TD> ddtrips, dptrips;
	for (int i = 0; i < this->refineMap.size(); i++){
		const int& i1 = this->refineMap[i].first;
		const int& i2 = this->refineMap[i].second;

		//quadratic self similarity
		double sim = aff.quadSelfSimDD(i1) + (i2 != -1 ? aff.quadSelfSimDD(i2) : 0) + (i2 != -1 ? 2.0*aff.simDD(i1, i2) : 0);
		if (sim > 1e-16){ddtrips.push_back(TD(i, i, sim));}

		//linear self similarity
		this->linaffdd(i) = aff.linSelfSim(i1) + (i2 != -1 ? aff.linSelfSim(i2) : 0);

		//data->data similarities
		for (int j = i+1; j < this->refineMap.size(); j++){
			const int& j1 = this->refineMap[j].first;
			const int& j2 = this->refineMap[j].second;
			double sim = aff.simDD(i1, j1) + (i2 != -1 ? aff.simDD(i2, j1) : 0) + (j2 != -1 ? aff.simDD(i1, j2) + (i2 != -1 ? aff.simDD(i2, j2) : 0) : 0);
			if (sim > 1e-16){
				ddtrips.push_back(TD(i, j, sim));
			}
		}
		//data->param similarities
		for (int j = 0; j < nOldPrm; j++){
			double sim = aff.simDP(i1, j) + (i2 != -1 ? aff.simDP(i2, j) : 0);
			if (sim > 1e-16){dptrips.push_back(TD(i, j, sim));}
		}
	}
	//param->param similarities
	for (int i = 0; i < nOldPrm; i++){
		this->affpp(i) = aff.selfSimPP(i, i);
	}
	//set the sparse matrices from triplets
	this->affdd.setFromTriplets(ddtrips.begin(), ddtrips.end());
	this->affdp.setFromTriplets(dptrips.begin(), dptrips.end());
	return;
}

template <class G>
double CoarseGraph<G>::linSelfSimDD(int i){
	return this->linaffdd(i);
}
template <class G>
double CoarseGraph<G>::quadSelfSimDD(int i){
	return this->affdd.coeff(i, i);
}
template <class G>
double CoarseGraph<G>::simDD(int i, int j){
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
double CoarseGraph<G>::selfSimPP(int i){
	return this->affpp(i);
}

template <class G>
double CoarseGraph<G>::simDP(int i, int j){
	return this->affdp.coeff(i, j);
}

template <class G>
std::vector<int> CoarseGraph<G>::getRefinedLabels(const std::vector<int>& lbls){
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
CoarseGraph<G>::getOldPrmLbls(){
	return this->oldPrmLbls;
}


//template <typename D, typename C, typename P>
//template <typename T> std::vector<int> KernDynMeans<D,C,P>::clusterSplit(std::vector<T>& data, std::vector<int> lbls){
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
//template <typename D, typename C, typename P>
//template <typename T> std::vector<int> KernDynMeans<D,C,P>::clusterMerge(std::vector<T>& data, std::vector<int> lbls){
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
