#ifndef __CLUSTER_IMPL_HPP

template<class D, class P>
Cluster<D, P>::Cluster(double lambda, double Q, double tau) : id(nextId++){
	this->age = 0;
	this->w = 0.0;
	this->lambda = lambda;
	this->Q = Q;
	this->tau = tau;
}

template<class D, class P>
double Cluster<D, P>::gamma() const{
	if (this->age == 0){ //age only =0 for brand new clusters
		return 0;
	} else {
		return 1.0/(1.0/this->w + this->tau*this->age);
	}
}


template<class D, class P>
void Cluster<D, P>::updatePrm(){
	this->prm.update(this->clusData.begin(), this->clusData.end(), this->gamma());
}

template<class D, class P>
std::vector<uint64_t> Cluster<D, P>::finalize(){
	std::vector<uint64_t> savedData;
	if(this->clusData.empty()){
		this->age++;
	} else {
		savedData = this->prm.updateOld(this->clusData.begin(), this->clusData.end(), this->gamma());
		this->w = this->gamma() + std::distance(this->clusData.begin(), this->clusData.end());
		this->age = 1;
	}
	this->clusData.clear();
	return savedData;
}

template<class D, class P>
std::vector<uint64_t> Cluster<D, P>::getAssignedIds() const{
	std::vector<uint64_t> asids;
	for (auto it = this->clusData.begin(); it != this->clusData.end(); ++it){
		asids.push_back(it->first);
	}
	return asids;
}

template<class D, class P>
void Cluster<D, P>::assignData(uint64_t did, D& d){
	if (this->clusData.find(did) != this->clusData.end()){
		throw DataAlreadyInClusterException(this->id, did);
	}
	this->clusData[did] = d;
}

template<class D, class P>
D Cluster<D, P>::deassignData(uint64_t did){
	if (this->clusData.find(did) == this->clusData.end()){
			throw DataNotInClusterException(this->id, did); 
	}
	D d = this->clusData[did];
	this->clusData.erase(did);
	return d;
}

template<class D, class P>
double Cluster<D, P>::distTo(const D& d) const{
	return this->prm.distTo(d, !this->clusData.empty());
}

template<class D, class P>
double Cluster<D, P>::cost() const{
	return this->clusData.empty() ? 0.0 : (this->age == 0 ? this->lambda : this->Q*this->age) 
		+ this->prm.cost(this->clusData.begin(), this->clusData.end(), this->gamma());
}

template<class D, class P>
bool Cluster<D, P>::isEmpty() const{
	return this->clusData.empty();
}

template<class D, class P>
bool Cluster<D, P>::isNew() const{
	return this->age == 0;
}

template<class D, class P>
P Cluster<D, P>::getPrm() const{
	return this->prm;
}

#define __CLUSTER_IMPL_HPP
#endif /* __CLUSTER_IMPL_HPP */
