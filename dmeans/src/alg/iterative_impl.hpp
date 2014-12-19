#ifndef __ITERATIVE_IMPL_HPP

template <class Model, bool monoCheck>
_Iterative<Model, monoCheck>::_Iterative(const Config& cfg){
	this->verbose = cfg.get("verbose", Config::OPTIONAL, false);
}

template <class Model, bool monoCheck>
double _Iterative<Model, monoCheck>::cluster(const std::vector<typename Model::Data>& obs, std::vector<Clus>& clus, const Model& model) const{
	if(obs.size() == 0){return 0.0;}
	//initial round of labelling data without deassigning it
	if(verbose){ std::cout << "Clustering step in iterative mode" << std::endl;}
	if(verbose){ std::cout << "Monotonicity checking: " << monoCheck << std::endl;}
	if(verbose){ std::cout << "Initial labelling..." << std::endl;}
	this->initialLabelling(obs, clus, model);
	//label/parameter update iteration
	bool labellingChanged = true;
	if (!monoCheck){
		while(labellingChanged){
			if(verbose){ std::cout << "Parameter update..." << std::endl;}
			this->parameterUpdate(clus, model);
			if(verbose){ std::cout << "Label update..." << std::endl;}
			labellingChanged = this->labelUpdate(clus, model);
		}
	} else {
		double obj = this->computeCost(clus, model);
		if(verbose){ std::cout << "Objective: " << obj << std::endl;}
		while(labellingChanged){
			double prevobj = obj;
			if(verbose){ std::cout << "Parameter update..." << std::endl;}
			this->parameterUpdate(clus, model);
			obj = this->computeCost(clus, model);
			if(verbose){ std::cout << "Objective: " << obj << std::endl;}
			if (obj > prevobj){
				throw MonotonicityViolationException(prevobj, obj, "parameterUpdate()");
			}
			prevobj = obj;
			if(verbose){ std::cout << "Label update..." << std::endl;}
			labellingChanged = this->labelUpdate(clus, model);
			obj = this->computeCost(clus, model);
			if(verbose){ std::cout << "Objective: " << obj << std::endl;}
			if (obj > prevobj){
				throw MonotonicityViolationException(prevobj, obj, "labelUpdate()");
			}
		}
	}
	if(verbose){ std::cout << "Done iterative clustering step" << std::endl;}
	return this->computeCost(clus, model);
}



template <class Model, bool monoCheck>
void _Iterative<Model, monoCheck>::initialLabelling(const std::vector< typename Model::Data>& obs, std::vector< Clus >& clus, const Model& model) const{
	std::vector<uint64_t> shuffs(obs.size());
	std::iota(shuffs.begin(), shuffs.end(), 0);
	std::shuffle(shuffs.begin(), shuffs.end(), RNG::get());
	for (uint64_t i = 0; i < shuffs.size(); i++){
		int id = shuffs[i];
		double minCost = std::numeric_limits<double>::infinity();
		uint64_t minInd = 0;
		for(uint64_t k = 0; k < clus.size(); k++){
			double d = model.compare(clus[k], obs[id]);
			if (d < minCost){
				minCost = d;
				minInd = k;
			}
		}
		if (model.exceedsNewClusterCost(obs[id], minCost)){
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
	std::shuffle(shuffs.begin(), shuffs.end(), RNG::get());
	//do reassignment
	for (uint64_t i = 0; i < shuffs.size(); i++){
		uint64_t lbl = lbls[shuffs[i]];
		uint64_t cid = lblToIdx[lbl];
		uint64_t id = ids[shuffs[i]];
		typename Model::Data obs = clus[cid].deassignData(id);
		bool deletedNew = false; //used to stop an infinite delete/create loop with labelchanged
		if (clus[cid].isNew() && clus[cid].isEmpty()){
			deletedNew = true;
			clus.erase(clus.begin()+cid);
			lblToIdx.erase(lbl);
			for(auto it = lblToIdx.begin(); it != lblToIdx.end(); ++it){
				if (it->second > cid){
					it->second--;
				}
			}
		}
		double minCost = std::numeric_limits<double>::infinity();
		uint64_t minInd = 0;
		for(uint64_t k = 0; k < clus.size(); k++){
			double d = model.compare(clus[k], obs);
			if (d < minCost){
				minCost = d;
				minInd = k;
			}
		}
		if (model.exceedsNewClusterCost(obs, minCost)){
			Clus newclus;
			newclus.assignData(id, obs);
			model.updatePrm(newclus);
			clus.push_back(newclus);
			if (!deletedNew){
				labellingChanged = true;
			}
		} else {
			if (cid != minInd){
				labellingChanged = true;
			}
			if (!clus[minInd].isNew() && clus[minInd].isEmpty()){
				clus[minInd].assignData(id, obs);
				model.updatePrm(clus[minInd]);
			} else {
				clus[minInd].assignData(id, obs);
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
