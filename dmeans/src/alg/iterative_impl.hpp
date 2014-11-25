#ifndef __ITERATIVE_IMPL_HPP

template <class Model, bool monoCheck>
_Iterative<Model, monoCheck>::_Iterative(const Config& cfg){
	this->cfg = cfg;
	this->verbose = cfg.get("verbose", Config::OPTIONAL, false);
}

template <class Model, bool monoCheck>
double _Iterative<Model, monoCheck>::cluster(const std::map<uint64_t, typename Model::Data>& obs, std::vector<Cluster<Model> >& clus, const Model& model) const{
	//initial round of labelling data without deassigning it
	this->initialLabelling(obs, clus);
	//label/parameter update iteration
	bool labellingChanged = true;
	if (!monoCheck){
		while(labellingChanged){
			this->parameterUpdate(clus);
			labellingChanged = this->labelUpdate(clus);
		}
	} else {
		double obj = this->computeCost(clus);
		while(labellingChanged){
			double prevobj = obj;
			std::cout << "prm update" << std::endl;
			this->parameterUpdate(clus);
			obj = this->computeCost(clus);
			std::cout << "Obj: " << obj << std::endl;
			if (obj > prevobj){
				throw MonotonicityViolationException(prevobj, obj, "parameterUpdate()");
			}
			prevobj = obj;
			std::cout << "lbl update" << std::endl;
			labellingChanged = this->labelUpdate(clus);
			obj = this->computeCost(clus);
			std::cout << "Obj: " << obj << std::endl;
			if (obj > prevobj){
				throw MonotonicityViolationException(prevobj, obj, "labelUpdate()");
			}
		}
	}
	return this->computeCost(clus);
}

template <class Model, bool monoCheck>
double _Iterative<Model, monoCheck>::computeCost(const std::vector< Cluster<Model> >& clus, const Model& model) const{
	double cost = 0;
	for(auto it = clus.begin(); it != clus.end(); ++it){
		cost += it->cost();
	}
	return cost;
}

template <class Model, bool monoCheck>
void _Iterative<Model, monoCheck>::initialLabelling(std::map<uint64_t, typename Model::Data>& obs, std::vector< Cluster<Model> >& clus, const Model& model){
	std::vector<uint64_t> ids;
	for (auto it = obs.begin(); it != obs.end(); ++it){
		ids.push_back(it->first);
	}
	std::random_shuffle(ids.begin(), ids.end());
	for (uint64_t i = 0; i < ids.size(); i++){
		double minCost = std::numeric_limits<double>::infinity();
		uint64_t minInd = -1;
		for(int k = 0; k < clus.size(); k++){
			double d = clus[k].compareTo(obs[ids[i]]);
			if (d < minCost){
				minCost = d;
				minInd = k;
			}
		}
		if (minCost > lambda){
			Cluster<Model> newclus(this->cfg);
			newclus.assignData(ids[i], obs[ids[i]]);
			newclus.updatePrm();
			clus.push_back(newclus);
		} else {
			clus[minInd].assignData(ids[i], obs[ids[i]]);
		}
	}
}

template <class Model, bool monoCheck>
bool _Iterative<Model, monoCheck>::labelUpdate(std::vector< Cluster<Model> >& clus, const Model& model){
	bool labellingChanged = false;
	//get the assignments across all clusters
	std::vector<uint64_t> ids, lbls;
	for (int k = 0; k < clus.size(); k++){
		std::vector<uint64_t> clusids = clus[k]getAssignedIds();
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
		typename Model::Data obs = clus[cl].deassignData(id);
		if (clus[cl].isNew() && clus[cl].isEmpty()){
			clus.erase(cl);
		}
		double minCost = std::numeric_limits<double>::infinity();
		uint64_t minInd = -1;
		for(int k = 0; k < clus.size(); k++){
			double d = clus[k].compareTo(obs);
			if (d < minCost){
				minCost = d;
				minInd = k;
			}
		}
		if (minCost > lambda){
			Cluster<Model> newclus(this->cfg);
			newclus.assignData(id, obs);
			newclus.updatePrm();
			clus.push_back(newclus);
			labellingChanged = true;
		} else {
			clus[minInd].assignData(id, obs);
			if (cl != minInd){
				labellingChanged = true;
			}
		}
	}
	return labellingChanged;
}

template <class Model, bool monoCheck>
void _Iterative<Model, monoCheck>::parameterUpdate(std::vector< Cluster<Model> >& clus, const Model& model){
	for (auto it = clus.begin(); it != clus.end(); ++it){
		it->updatePrm();
	}
}

#define __ITERATIVE_IMPL_HPP
#endif 
