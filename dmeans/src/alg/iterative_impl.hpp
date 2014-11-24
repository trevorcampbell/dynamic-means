#ifndef __ITERATIVE_IMPL_HPP
template <class D, class P, bool M>
void _Iterative<D, P, M>::reset(){
	this->clusters.clear();
}

template <class D, class P, bool M>
void _Iterative<D, P, M>::finalize(){
	for(uint64_t k = 0; k < this->clusters.size(); k++){
		this->clusters[k].finalize();
	}
}

template <class D, class P, bool M>
void _Iterative<D, P, M>::cluster(std::map<uint64_t, D>& obs, std::map<uint64_t, Cluster<D, P> >& clus, double lambda, double Q, double tau, bool verbose){
	this->timer.start();
	std::vector<Cluster<D, P> > savedInitialClusters = this->clusters;
	Results<P> bestResults;
	std::vector<Cluster<D, P> > bestClusters;
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
			bestClusters = this->clusters;
			bestResults.obj = obj;
		}
		this->clusters = savedInitialClusters;
	}
	this->clusters = bestClusters;
	this->finalize();
	bestResults.tTaken = this->timer.elapsed_ms();
	return bestResults;
}

template <class D, class P, bool M>
double _Iterative<D, P, M>::computeCost(){
	double cost = 0;
	for(uint64_t k = 0; k < this->clusters.size(); k++){
		cost += this->clusters[k].cost();
	}
	return cost;
}

template <class D, class P, bool M>
void _Iterative<D, P, M>::initialLabelling(std::map<uint64_t, D>& obs){
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

template <class D, class P, bool M>
bool _Iterative<D, P, M>::labelUpdate(){
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

template <class D, class P, bool M>
void _Iterative<D, P, M>::parameterUpdate(){
	for (uint64_t k = 0; k < this->clusters.size(); k++){
		this->clusters[k].updatePrm();
	}
}

#define __ITERATIVE_IMPL_HPP
#endif 
