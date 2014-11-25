#ifndef __DMEANS_IMPL_HPP
template <class Model, template<typename> class Alg>
DMeans<Model, Alg>::DMeans(Config cfg){
	this->cfg = cfg;
	this->verbose = cfg.get("verbose", Config::OPTIONAL, false);
	this->nextLabel = 0;
	uint64_t seed = cfg.get("seed", Config::OPTIONAL, -1);
	if (seed < 0){
		std::srand(this->timer.now_ms());
	} else {
		std::srand(seed);
	}
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
		if (it->isPermanentlyDead()){
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
Results<Model> DMeans<Model, Alg>::cluster(std::vector<typename Model::Data>& obs, uint64_t nRestarts){
	this->timer.start();//start the timer
	//these store the best cost/clustering over nRestarts restarts
	double minCost = std::numeric_limits<double>::infinity();
	std::vector<Cluster<Model> > minClusters;
	//convert the data to a map with unique indices for insertion into clusters
	std::map<uint64_t, typename Model::Data> obsMap;
	for(uint64_t i = 0; i < obs.size(); i++){
		obsMap[i] = obs[i];
	}
	Alg<Model> alg(this->cfg);
	for (uint64_t k = 0; k < nRestarts; k++){
		//cluster at this step
		double cost = alg.cluster(obsMap, this->clusters, this->verbose);
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
	bestResults.tTaken = this->timer.elapsed_ms();
	this->finalize();
	return bestResults;
}

template <class Model, template<typename> class Alg>
Results<Model> DMeans<Model, Alg>::getResults() const{
	Results<Model> r;
	for(auto it = this->clusters.begin(); it != this->clusters.end(); ++it){
		r.prms[it->id] = it->getPrm();
		std::vector<uint64_t> dids = it->getAssignedIds();
		for(uint64_t i = 0; i < dids.size(); i++){
			r.lbls[dids[i]] = it->id;
		}
	}
	return r;
}


#define __DMEANS_IMPL_HPP
#endif 
