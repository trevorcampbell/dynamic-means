#ifndef __CLUSTER_IMPL_HPP

template<class P, class D>
Cluster<P, D>::Cluster() : id(nextId++){
	this->age = 0;
	this->w = 0.0;
}

template<class P, class D>
double Cluster<P, D>::gamma(double tau){
	if (this->age == 0){ //age only =0 for brand new clusters
		return 0;
	} else {
		return 1.0/(1.0/this->w + tau*this->age);
	}
}


template<class P, class D>
void Cluster<P, D>::updatePrm(){
	this->prm.update(this->clusData.begin(), this->clusData.end(), this->gamma());
}

template<class P, class D>
std::vector<uint64_t> Cluster<P, D>::finalize(){
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

template<class P, class D>
std::vector<uint64_t> Cluster<P, D>::getAssignedIDs() const{
	std::vector<uint64_t> asids;
	for (auto it = this->clusData.begin(), it != this->clusData.end(); ++it){
		asids.push_back(it->first);
	}
	return asids;
}

template<class P, class D>
void Cluster<P, D>::assignData(Data<D>& d){
	if (this->clusData.find(d.id) != this->clusData.end()){
		throw DataAlreadyInClusterException(this->id, d.id);
	}
	this->clusData[d.id] = d;
}

template<class P, class D>
Data<D> Cluster<P, D>::deassignData(uint64_t did){
	if (this->clusData.find(did) == this->clusData.end()){
			throw DataNotInClusterException(this->id, did); 
	}
	Data d = this->clusData[did];
	this->clusData.erase(did);
	return d;
}

template<class P, class D>
double Cluster<P, D>::distTo(const Data<D>& d) const{
	return this->prm.distTo(d.d, !this->clusData.empty());
}

template<class P, class D>
double Cluster<P, D>::cost(double lambda, double Q) const{
	return this->clusData.empty() ? 0.0 : (this->age == 0 ? lambda : Q*this->age) + this->prm.cost(this->clusData.begin(), this->clusData.end(), this->gamma());
}

template<class P, class D>
bool Cluster<P, D>::isEmpty() const{
	return this->clusData.empty();
}

template<class P, class D>
bool Cluster<P, D>::isNew() const{
	return this->age == 0;
}

template<class P, class D>
P Cluster<P, D>::getPrm() const{
	return this->prm;
}

#define __CLUSTER_IMPL_HPP
#endif /* __CLUSTER_IMPL_HPP */
