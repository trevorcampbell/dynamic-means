#ifndef __DMEANS_IMPL_HPP
template <class Model, template<typename> class Alg>
DMeans<Model, Alg>::DMeans(Config cfg) : model(cfg){
	this->cfg = cfg;
	this->verbose = cfg.get("verbose", Config::Type::OPTIONAL, false);
	this->nRestarts = this->cfg.get("nRestarts", Config::Type::OPTIONAL, 1);
	this->nextLabel = 0;
	int s = cfg.get("seed", Config::Type::OPTIONAL, -1);
	if (s >= 0){ RNG::get().seed(s); }
}

template <class Model, template<typename> class Alg>
void DMeans<Model, Alg>::reset(){
	this->clusters.clear();
	this->nextLabel = 0;
}

template <class Model, template<typename> class Alg>
void DMeans<Model, Alg>::finalize(){
	auto it = this->clusters.begin();
	while(it != this->clusters.end()){
		it->finalize();
		if(model.isClusterDead(it->getAge())){
			it = this->clusters.erase(it);
		} else {
			++it;
		}
	}
}

template <class Model, template<typename> class Alg>
void DMeans<Model, Alg>::restart(){
	auto it = this->clusters.begin();
	while(it != this->clusters.end()){
		if (it->isNew()){ //remove new clusters entirely
			it = this->clusters.erase(it);
		} else { //erase data from old clusters
			it->clearData();
			++it;
		}
	}
}

template <class Model, template<typename> class Alg>
void DMeans<Model, Alg>::labelNewClusters(){
	for(auto it = this->clusters.begin(); it != this->clusters.end(); ++it){
		if (it->isNew()){
			it->setID(this->nextLabel++);
		}
	}
}


template <class Model, template<typename> class Alg>
Results<Model> DMeans<Model, Alg>::cluster(const std::vector<typename Model::Data>& obs){
	this->timer.start();//start the timer
	//these store the best cost/clustering over nRestarts restarts
	double minCost = std::numeric_limits<double>::infinity();
	std::vector<Cluster<typename Model::Data, typename Model::Parameter> > minClusters;
	Alg<Model> alg(this->cfg);
	for (uint64_t k = 0; k < this->nRestarts; k++){
		//cluster at this step
		double cost = alg.cluster(obs, this->clusters, this->model);
		if (cost < minCost){
			//the objective is the best so far, save the results
			minCost = cost;
			minClusters = this->clusters;
		}
		//reset the clusters to their initial state for this step
		this->restart();
	}
	//restore the best clustering
	this->clusters = minClusters;
	//add labels to all new clusters
	this->labelNewClusters();
	//output the results
	Results<Model> bestResults = this->getResults();
	bestResults.cost = minCost;
	bestResults.tTaken = this->timer.elapsed_s();
	this->finalize();
	return bestResults;
}

template <class Model, template<typename> class Alg>
Results<Model> DMeans<Model, Alg>::getResults() const{
	Results<Model> r;
	for(auto it = this->clusters.begin(); it != this->clusters.end(); ++it){
		r.prms[it->getID()] = it->getPrm();
		std::vector<uint64_t> dids = it->getAssignedIds();
		for(uint64_t i = 0; i < dids.size(); i++){
			r.lbls[dids[i]] = it->getID();
		}
	}
	return r;
}


#define __DMEANS_IMPL_HPP
#endif 
