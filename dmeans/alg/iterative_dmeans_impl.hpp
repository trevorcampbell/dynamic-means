#ifndef __ITERATIVE_DMEANS_IMPL_HPP
template <class D, class P>
IterativeDMeans<D, P>::IterativeDMeans(double lambda, double Q, double tau, bool verbose, int seed){
	this->lambda = lambda;
	this->Q = Q;
	this->tau = tau;
	this->verbose = verbose;
	if (seed < 0){
		std::srand(this->timer.now_ms());
	} else {
		std::srand(seed);
	}
}

template <class D, class P>
void IterativeDMeans<D, P>::reset(){
	this->clusters.clear();
}

template <class D, class P>
Results<P> IterativeDMeans<D, P>::cluster(std::map<uint64_t, D>& obs, uint64_t nRestarts, bool checkCosts){
	this->timer.start();
	Results<P> bestResults;
	bestResults.obj = std::numeric_limits<double>::infinity();
	for (uint64_t restart = 0; restart < nRestarts; restart++){
		//initial round of labelling data without deassigning it
		this->initialLabelling(obs);
		//label/parameter update iteration
		bool labellingChanged = true;
		if (!checkCosts){
			while(labellingChanged){
				this->parameterUpdate();
				labellingChanged = this->labelUpdate();
			}
		} else {
			double obj = this->computeCost();
			while(labellingChanged){
				double prevobj = obj;
				this->parameterUpdate();
				obj = this->computeCost();
				if (obj > prevobj){
					throw MonotonicityViolationException(prevobj, obj, "parameterUpdate()");
				}
				prevobj = obj;
				labellingChanged = this->labelUpdate();
				obj = this->computeCost();
				if (obj > prevobj){
					throw MonotonicityViolationException(prevobj, obj, "labelUpdate()");
				}
			}
		}
		double obj = this->computeCost();
		if (obj < bestResults.obj){
			//the objective is the best so far, save the results
			bestResults = this->computeResults();
			bestResults.obj = obj;
		}
	}
	bestResults.tTaken = this->timer.elapsed_ms();
	return bestResults;
}

template <class D, class P>
Results<P> IterativeDMeans<D, P>::computeResults(){
	Results<P> r;
	for(uint64_t k = 0; k < this->clusters.size(); k++){
		r.prms[this->clusters[k].id] = this->clusters[k].getPrm();
		std::vector<uint64_t> dids = this->clusters[k].getAssignedIds();
		for(uint64_t i = 0; i < dids.size(); i++){
			r.lbls[dids[i]] = this->clusters[k].id;
		}
	}
	return r;
}

template <class D, class P>
double IterativeDMeans<D, P>::computeCost(){
	double cost = 0;
	for(uint64_t k = 0; k < this->clusters.size(); k++){
		cost += this->clusters[k].cost();
	}
	return cost;
}

template <class D, class P>
void IterativeDMeans<D, P>::initialLabelling(std::map<uint64_t, D>& obs){
	std::vector<uint64_t> ids;
	for (auto it = obs.begin(); it != obs.end(); ++it){
		ids.push_back(it->first);
	}
	std::random_shuffle(ids.begin(), ids.end());
	for (uint64_t i = 0; i < ids.size(); i++){
		double minCost = std::numeric_limits<double>::infinity();
		uint64_t minInd = -1;
		for(uint64_t k = 0; k < this->clusters.size(); k++){
			double d = this->clusters[k].distTo(obs[ids[i]]);
			if (d < minCost){
				minCost = d;
				minInd = k;
			}
		}
		if (minCost > this->lambda){
			this->clusters.push_back(Cluster<D, P>(this->lambda, this->Q, this->tau));
			this->clusters.back().assignData(ids[i], obs[ids[i]]);
		} else {
			this->clusters[minInd].assignData(ids[i], obs[ids[i]]);
		}
	}
}

template <class D, class P>
bool IterativeDMeans<D, P>::labelUpdate(){
	bool labellingChanged = false;
	//get the assignments across all clusters
	std::vector<uint64_t> ids, lbls;
	for (uint64_t k = 0; k < this->clusters.size(); k++){
		std::vector<uint64_t> clusids = this->clusters[k].getAssignedIds();
		ids.insert(ids.end(), clusids.begin(), clusids.end());
		lbls.insert(lbls.end(), clusids.size(), k);
	}
	std::vector<uint64_t> shuffs(ids.size());
	std::iota(shuffs.begin(), shuffs.end(), 0);
	std::random_shuffle(shuffs.begin(), shuffs.end());
	//do reassignment
	for (uint64_t i = 0; i < shuffs.size(); i++){
		uint64_t cl = lbls[shuffs[i]];
		uint64_t id = ids[shuffs[i]];
		D obs = this->clusters[cl].deassignData(id);
		if (this->clusters[cl].isNew() && this->clusters[cl].isEmpty()){
			this->clusters.erase(this->clusters.begin()+cl);
		}
		double minCost = std::numeric_limits<double>::infinity();
		uint64_t minInd = -1;
		for(uint64_t k = 0; k < this->clusters.size(); k++){
			double d = this->clusters[k].distTo(obs);
			if (d < minCost){
				minCost = d;
				minInd = k;
			}
		}
		if (minCost > this->lambda){
			this->clusters.push_back(Cluster<D, P>(this->lambda, this->Q, this->tau));
			this->clusters.back().assignData(id, obs);
			this->clusters.back().updatePrm();
			labellingChanged = true;
		} else {
			this->clusters[minInd].assignData(id, obs);
			labellingChanged = true;
		}
	}
	return labellingChanged;
}

template <class D, class P>
void IterativeDMeans<D, P>::parameterUpdate(){
	for (uint64_t k = 0; k < this->clusters.size(); k++){
		this->clusters[k].updatePrm();
	}
}

#define __ITERATIVE_DMEANS_IMPL_HPP
#endif /* __ITERATIVE_DMEANS_IMPL_HPP */
