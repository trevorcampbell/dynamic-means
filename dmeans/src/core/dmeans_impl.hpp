#ifndef __ITERATIVE_DMEANS_IMPL_HPP
template <class D, class P, class A>
DMeans<D, P, A>::DMeans(double lambda, double Q, double tau, bool verbose, int seed){
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

template <class D, class P, class A>
void DMeans<D, P, A>::reset(){
	this->clusters.clear();
}

template <class D, class P, class A>
void DMeans<D, P, A>::finalize(){
	for(auto it = this->clusters.begin(); it != this->clusters.end(); ++it){
		it->second.finalize(this->tau);
	}
}


template <class D, class P, class A>
void DMeans<D, P, A>::restart(){
	for(auto it = this->clusters.begin(); it != this->clusters.end(); ++it){
		it->second.clearData();
	}
}

template <class D, class P, class A>
Results<P> DMeans<D, P, A>::cluster(std::map<uint64_t, D>& obs, uint64_t nRestarts){
	this->timer.start();
	Results<P> bestResults;
	bestResults.cost = std::numeric_limits<double>::infinity();
	for (uint64_t k = 0; k < nRestarts; k++){
		//cluster at this step
		A::cluster(obs, this->clusters, this->lambda, this->Q, this->tau);
		//compute the resulting objective
		Results<P> res = this->computeResults();
		if (res.cost < bestResults.cost){
			//the objective is the best so far, save the results
			bestResults = res;
		}
		//reset the clusters to their initial state for this step
		this->restart();
	}
	//restore the best clustering
	for (auto it = obs.begin(); it != obs.end(); ++it){
		uint64_t lbl = bestResults.lbls[it->first];
		this->clusters[lbl].assignData(it->first, it->second);
	}
	this->finalize();
	bestResults.tTaken = this->timer.elapsed_ms();
	return bestResults;
}

template <class D, class P, class A>
Results<P> DMeans<D, P, A>::computeResults() const{
	Results<P> r;
	r.cost = 0;
	for(auto it = this->clusters.begin(); it != this->clusters.end(); ++it){
		r.prms[it->first] = it->second.getPrm();
		std::vector<uint64_t> dids = it->second.getAssignedIds();
		for(uint64_t i = 0; i < dids.size(); i++){
			r.lbls[dids[i]] = it->first;
		}
		r.cost += it->second.cost(this->lambda, this->Q);
	}
	return r;
}


#define __ITERATIVE_DMEANS_IMPL_HPP
#endif /* __ITERATIVE_DMEANS_IMPL_HPP */
