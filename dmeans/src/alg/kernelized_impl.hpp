#ifndef __KERNELIZED_IMPL_HPP

template <class Model, bool monoCheck>
_Kernelized<Model, monoCheck>::_Kernelized(const Config& cfg){
	this->verbose = cfg.get("verbose", Config::OPTIONAL, false);
}

template <class Model, bool monoCheck>
double _Kernelized<Model, monoCheck>::cluster(const std::vector<typename Model::Data>& obs, std::vector<Clus>& clus, const Model& model) const{
	if(obs.size() == 0){return 0.0;}
	//initial round of labelling data without deassigning it
	if(verbose){ std::cout << "Clustering step in kernelized mode" << std::endl;}
	if(verbose){ std::cout << "Monotonicity checking: " << monoCheck << std::endl;}
	if(verbose){ std::cout << "Initial labelling..." << std::endl;}
	this->initialLabelling(obs, clus, model);
	//label/parameter update iteration
	bool labellingChanged = true;
	if (!monoCheck){
		while(labellingChanged){
			if(verbose){ std::cout << "Label update..." << std::endl;}
			labellingChanged = this->labelUpdate(clus, model);
		}
	} else {
		double obj = this->computeCost(clus, model);
		if(verbose){ std::cout << "Objective: " << obj << std::endl;}
		while(labellingChanged){
			double prevobj = obj;
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
	this->parameterUpdate(clus, model);
	return this->computeCost(clus, model);
}



template <class Model, bool monoCheck>
void _Kernelized<Model, monoCheck>::initialLabelling(const std::vector< typename Model::Data>& obs, std::vector< Clus >& clus, const Model& model) const{
	std::vector<uint64_t> toAsgn(obs.size(), 0);
	std::iota(toAsgn.begin(), toAsgn.end(), 0);
	std::uniform_int_distribution<> uid(0, obs.size()-1);
	int id = uid(RNG::get());
	Clus newclus;
	model.kernelAdd(newclus, obs[id]);
	newclus.assignData(id, obs[id]);
	clus.push_back(newclus);
	toAsgn.erase(toAsgn.begin()+id);
	while (!toAsgn.empty()){
		double minCost = std::numeric_limits<double>::infinity();
		uint64_t minK = 0;
		uint64_t minInd = 0;
		for (uint64_t i = 0; i < toAsgn.size(); i++){
			int id = toAsgn[i];
			for(uint64_t k = 0; k < clus.size(); k++){
				double d = model.compare(clus[k], obs[id]);
				if (d < minCost){
					minCost = d;
					minK = k;
					minInd = i;
				}
			}
		}
		int id = toAsgn[minInd];
		if (model.exceedsNewClusterCost(obs[id], minCost)){
			Clus newclus;
			model.kernelAdd(newclus, obs[id]);
			newclus.assignData(id, obs[id]);
			clus.push_back(newclus);
		} else {
			model.kernelAdd(clus[minK], obs[id]);
			clus[minK].assignData(id, obs[id]);
		}
		toAsgn.erase(toAsgn.begin()+minInd);
	}
}

template <class Model, bool monoCheck>
bool _Kernelized<Model, monoCheck>::labelUpdate(std::vector< Clus >& clus, const Model& model) const{
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
		model.kernelSubtract(clus[cid], obs);
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
			model.kernelAdd(newclus, obs);
			newclus.assignData(id, obs);
			clus.push_back(newclus);
			if (!deletedNew){
				labellingChanged = true;
			}
		} else {
			if (cid != minInd){
				labellingChanged = true;
			}
			model.kernelAdd(clus[minInd], obs);
			clus[minInd].assignData(id, obs);
		}
	}
	return labellingChanged;
}

template <class Model, bool monoCheck>
void _Kernelized<Model, monoCheck>::parameterUpdate(std::vector< Clus >& clus, const Model& model) const{
	for (auto it = clus.begin(); it != clus.end(); ++it){
		model.updatePrm(*it);
	}
}

template <class Model, bool monoCheck>
double _Kernelized<Model, monoCheck>::computeCost(const std::vector< Clus >& clus, const Model& model) const{
	double cost = 0;
	for(auto it = clus.begin(); it != clus.end(); ++it){
		cost += model.clusterCost(*it);
	}
	return cost;
}

#define __KERNELIZED_IMPL_HPP
#endif 
