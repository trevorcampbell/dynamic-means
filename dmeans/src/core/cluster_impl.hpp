#ifndef __CLUSTER_IMPL_HPP

template<class D, class P>
Cluster<D, P>::Cluster(){
	this->age = 0;
	this->idset = false;
	this->id = 0;
}


template<class D, class P>
Cluster<D, P>::Cluster(const Cluster<D, P>& rhs) {
	this->age = rhs.age;
	this->prm = rhs.prm;
	this->oldprm = rhs.oldprm;
	this->clusData = rhs.clusData;
	this->id = rhs.id;
	this->idset = rhs.idset;
}


template<class D, class P>
Cluster<D, P>& Cluster<D, P>::operator=(const Cluster<D, P>& rhs) {
	if (this != &rhs){
		this->age = rhs.age;
		this->prm = rhs.prm;
		this->oldprm = rhs.oldprm;
		this->clusData = rhs.clusData;
		this->id = rhs.id;
		this->idset = rhs.idset;
	}
	return *this;
}

template<class D, class P>
void Cluster<D, P>::finalize(){
	if(this->isEmpty()){
		this->age++;
	} else {
		this->oldprm = this->prm;
		this->age = 1;
	}
	this->clusData.clear();
}

template<class D, class P>
uint64_t Cluster<D, P>::getAge() const{
	return this->age;
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
void Cluster<D, P>::assignData(uint64_t did, const D& d){
	if (this->clusData.find(did) != this->clusData.end()){
		throw DataAlreadyInClusterException(did);
	}
	this->clusData[did] = d;
}

template<class D, class P>
D Cluster<D, P>::deassignData(uint64_t did){
	if (this->clusData.find(did) == this->clusData.end()){
		throw DataNotInClusterException(did); 
	}
	D d = this->clusData[did];
	this->clusData.erase(did);
	return d;
}

template<class D, class P>
void Cluster<D, P>::clearData(){
	this->clusData.clear();
}

template<class D, class P>
void Cluster<D, P>::setID(uint64_t id){
	if (this->idset){
		throw IDAlreadySetException(this->id, id);
	}
	this->id = id;
	this->idset = true;
}

template<class D, class P>
uint64_t Cluster<D, P>::getID() const{
	if (!this->idset){
		throw IDUnsetException();
	}
	return this->id;
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
P& Cluster<D, P>::getPrmRef(){
	return this->prm;
}

template<class D, class P>
P& Cluster<D, P>::getOldPrmRef(){
	return this->oldprm;
}

template<class D, class P>
const P& Cluster<D, P>::getPrm() const{
	return this->prm;
}

template<class D, class P>
const P& Cluster<D, P>::getOldPrm() const{
	return this->oldprm;
}

template<class D, class P>
typename std::map<uint64_t, D>::const_iterator Cluster<D, P>::data_cbegin() const{
	return this->clusData.cbegin();
}

template<class D, class P>
typename std::map<uint64_t, D>::const_iterator Cluster<D, P>::data_cend() const{
	return this->clusData.cend();
}

#define __CLUSTER_IMPL_HPP
#endif /* __CLUSTER_IMPL_HPP */
