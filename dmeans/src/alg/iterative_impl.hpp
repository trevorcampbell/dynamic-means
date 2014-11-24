#ifndef __ITERATIVE_IMPL_HPP
template <class D, class P, bool M>
void _Iterative<D, P, M>::cluster(std::map<uint64_t, D>& obs, std::map<uint64_t, Cluster<D, P> >& clus, double lambda, double Q, double tau, bool verbose){
	//initial round of labelling data without deassigning it
	this->initialLabelling(obs, clus, lambda);
	//label/parameter update iteration
	bool labellingChanged = true;
	if (!M){
		while(labellingChanged){
			this->parameterUpdate(clus);
			labellingChanged = this->labelUpdate(clus, lambda);
		}
	} else {
		double obj = this->computeCost(clus, lambda, Q);
		while(labellingChanged){
			double prevobj = obj;
			this->parameterUpdate(clus);
			obj = this->computeCost(clus, lambda, Q);
			if (obj > prevobj){
				throw MonotonicityViolationException(prevobj, obj, "parameterUpdate()");
			}
			prevobj = obj;
			labellingChanged = this->labelUpdate(clus, lambda);
			obj = this->computeCost(clus, lambda, Q);
			if (obj > prevobj){
				throw MonotonicityViolationException(prevobj, obj, "labelUpdate()");
			}
		}
	}
}

template <class D, class P, bool M>
double _Iterative<D, P, M>::computeCost(std::map<uint64_t, Cluster<D, P> >& clus, double lambda, double Q){
	double cost = 0;
	for(auto it = clus.begin(); it != clus.end(); ++it){
		cost += it->second.cost(lambda, Q);
	}
	return cost;
}

template <class D, class P, bool M>
void _Iterative<D, P, M>::initialLabelling(std::map<uint64_t, D>& obs, std::map<uint64_t, Cluster<D, P> >& clus, double lambda){
	std::vector<uint64_t> ids;
	for (auto it = obs.begin(); it != obs.end(); ++it){
		ids.push_back(it->first);
	}
	std::random_shuffle(ids.begin(), ids.end());
	for (uint64_t i = 0; i < ids.size(); i++){
		double minCost = std::numeric_limits<double>::infinity();
		uint64_t minInd = -1;
		for(auto it = clus.begin(); it != clus.end(); ++it){
			double d = it->second.distTo(obs[ids[i]]);
			if (d < minCost){
				minCost = d;
				minInd = it->first;
			}
		}
		if (minCost > lambda){
			Cluster<D, P> newclus;
			newclus.assignData(ids[i], obs[ids[i]]);
			newclus.updatePrm();
			clus[newclus.id()] = newclus;
		} else {
			clus[minInd].assignData(ids[i], obs[ids[i]]);
		}
	}
}

template <class D, class P, bool M>
bool _Iterative<D, P, M>::labelUpdate(std::map<uint64_t, Cluster<D, P> >& clus, double lambda){
	bool labellingChanged = false;
	//get the assignments across all clusters
	std::vector<uint64_t> ids, lbls;
	for (auto it = clus.begin(); it != clus.end(); ++it){
		std::vector<uint64_t> clusids = it->second.getAssignedIds();
		ids.insert(ids.end(), clusids.begin(), clusids.end());
		lbls.insert(lbls.end(), clusids.size(), it->first);
	}
	std::vector<uint64_t> shuffs(ids.size());
	std::iota(shuffs.begin(), shuffs.end(), 0);
	std::random_shuffle(shuffs.begin(), shuffs.end());
	//do reassignment
	for (uint64_t i = 0; i < shuffs.size(); i++){
		uint64_t cl = lbls[shuffs[i]];
		uint64_t id = ids[shuffs[i]];
		D obs = clus[cl].deassignData(id);
		if (clus[cl].isNew() && clus[cl].isEmpty()){
			clus.erase(cl);
		}
		double minCost = std::numeric_limits<double>::infinity();
		uint64_t minInd = -1;
		for(auto it = clus.begin(); it != clus.end(); ++it){
			double d = it->second.distTo(obs);
			if (d < minCost){
				minCost = d;
				minInd = it->first;
			}
		}
		if (minCost > lambda){
			Cluster<D, P> newclus;
			newclus.assignData(id, obs);
			newclus.updatePrm();
			clus[newclus.id()] = newclus;
			labellingChanged = true;
		} else {
			clus[minInd].assignData(id, obs);
			labellingChanged = true;
		}
	}
	return labellingChanged;
}

template <class D, class P, bool M>
void _Iterative<D, P, M>::parameterUpdate(std::map<uint64_t, Cluster<D, P> >& clus){
	for (auto it = clus.begin(); it != clus.end(); ++it){
		it->second.updatePrm();
	}
}

#define __ITERATIVE_IMPL_HPP
#endif 
