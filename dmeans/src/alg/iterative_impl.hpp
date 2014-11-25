#ifndef __ITERATIVE_IMPL_HPP

template <class Model, bool monoCheck>
_Iterative<Model, monoCheck>::_Iterative(const Config& cfg){
	this->cfg = cfg;
	this->verbose = cfg.get("verbose", Config::OPTIONAL, false);
}

template <class Model, bool monoCheck>
double _Iterative<Model, monoCheck>::cluster(const std::vector<typename Model::Data>& obs, std::vector<Clus>& clus, const Model& model) const{
	std::cout << "Clustering " << obs.size() << " datapoints" << std::endl;
	std::cout << "Initial # clusters: " << clus.size() << std::endl;
	//initial round of labelling data without deassigning it
	this->initialLabelling(obs, clus, model);
	std::cout  << "After initial labelling # clusters: " << clus.size() << std::endl;
	//label/parameter update iteration
	bool labellingChanged = true;
	if (!monoCheck){
		while(labellingChanged){
			this->parameterUpdate(clus, model);
			labellingChanged = this->labelUpdate(clus, model);
		}
	} else {
		double obj = this->computeCost(clus, model);
		std::cout << "Obj: " << obj << std::endl;
		while(labellingChanged){
			double prevobj = obj;
			std::cout << "prmupdate" << std::endl;
			this->parameterUpdate(clus, model);
			obj = this->computeCost(clus, model);
			std::cout << "Obj: " << obj << std::endl;
			if (obj > prevobj){
				throw MonotonicityViolationException(prevobj, obj, "parameterUpdate()");
			}
			prevobj = obj;
			std::cout << "lblupdate" << std::endl;
			labellingChanged = this->labelUpdate(clus, model);
			obj = this->computeCost(clus, model);
			std::cout << "Obj: " << obj << std::endl;
			if (obj > prevobj){
				throw MonotonicityViolationException(prevobj, obj, "labelUpdate()");
			}
		}
	}
	return this->computeCost(clus, model);
}



template <class Model, bool monoCheck>
void _Iterative<Model, monoCheck>::initialLabelling(const std::vector< typename Model::Data>& obs, std::vector< Clus >& clus, const Model& model) const{
	std::vector<uint64_t> shuffs(obs.size());
	std::iota(shuffs.begin(), shuffs.end(), 0);
	std::random_shuffle(shuffs.begin(), shuffs.end());
	for (uint64_t i = 0; i < shuffs.size(); i++){
		int id = shuffs[i];
		double minCost = std::numeric_limits<double>::infinity();
		uint64_t minInd = -1;
		for(uint64_t k = 0; k < clus.size(); k++){
			double d = model.compare(clus[k], obs[id]);
			if (d < minCost){
				minCost = d;
				minInd = k;
			}
		}
		if (model.exceedsNewClusterCost(obs[id], minCost)){
			std::cout << "Creating cluster in initial labelling step " << std::endl;
			Clus newclus;
			newclus.assignData(id, obs[id]);
			model.updatePrm(newclus);
			clus.push_back(newclus);
		} else {
			clus[minInd].assignData(id, obs[id]);
		}
	}
}

template <class Model, bool monoCheck>
bool _Iterative<Model, monoCheck>::labelUpdate(std::vector< Clus >& clus, const Model& model) const{
	bool labellingChanged = false;
	//first create a map of lbl to cluster index to make deletion easier
	std::map<uint64_t, uint64_t> lblToIdx;
	for (uint64_t k = 0; k < clus.size(); k++){
		lblToIdx[k] = k;
	}

	//get the assignments across all clusters
	std::vector<uint64_t> ids, lbls;
	for (uint64_t k = 0; k < clus.size(); k++){
		std::vector<uint64_t> clusids = clus[k].getAssignedIds();
		ids.insert(ids.end(), clusids.begin(), clusids.end());
		lbls.insert(lbls.end(), clusids.size(), k);
	}
	std::vector<uint64_t> shuffs(ids.size());
	std::iota(shuffs.begin(), shuffs.end(), 0);
	std::random_shuffle(shuffs.begin(), shuffs.end());
	//do reassignment
	for (uint64_t i = 0; i < shuffs.size(); i++){
		uint64_t lbl = lbls[shuffs[i]];
		uint64_t cid = lblToIdx[lbl];
		uint64_t id = ids[shuffs[i]];
		typename Model::Data obs = clus[cid].deassignData(id);
		if (clus[cid].isNew() && clus[cid].isEmpty()){
			std::cout << "Deleting cluster at idx " << cid << " in initial labelling step " << std::endl;
			clus.erase(clus.begin()+cid);
			lblToIdx.erase(lbl);
			for(auto it = lblToIdx.begin(); it != lblToIdx.end(); ++it){
				if (it->second > cid){
					it->second--;
				}
			}
		}
		double minCost = std::numeric_limits<double>::infinity();
		uint64_t minInd = -1;
		for(uint64_t k = 0; k < clus.size(); k++){
			double d = model.compare(clus[k], obs);
			if (d < minCost){
				minCost = d;
				minInd = k;
			}
		}
		if (model.exceedsNewClusterCost(obs, minCost)){
			std::cout << "Creating cluster in initial labelling step " << std::endl;
			Clus newclus;
			newclus.assignData(id, obs);
			model.updatePrm(newclus);
			clus.push_back(newclus);
			labellingChanged = true;
		} else {
			clus[minInd].assignData(id, obs);
			if (cid != minInd){
				labellingChanged = true;
			}
		}
	}
	return labellingChanged;
}

template <class Model, bool monoCheck>
void _Iterative<Model, monoCheck>::parameterUpdate(std::vector< Clus >& clus, const Model& model) const{
	for (auto it = clus.begin(); it != clus.end(); ++it){
		model.updatePrm(*it);
	}
}

template <class Model, bool monoCheck>
double _Iterative<Model, monoCheck>::computeCost(const std::vector< Clus >& clus, const Model& model) const{
	double cost = 0;
	for(auto it = clus.begin(); it != clus.end(); ++it){
		cost += model.clusterCost(*it);
	}
	return cost;
}

#define __ITERATIVE_IMPL_HPP
#endif 
