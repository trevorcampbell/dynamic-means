#ifndef __CLUSTER_IMPL_HPP

template<class Model>
Cluster<Model>::Cluster(){
	this->age = 0;
	this->w = 0.0;
	this->gamma = 0.0;
	this->id = -1;
}


template<class Model>
Cluster<Model>::Cluster(const Cluster<Model>& rhs) {
	this->age = rhs.age;
	this->w = rhs.w;
	this->gamma = rhs.gamma;
	this->prm = rhs.prm;
	this->clusData = rhs.clusData;
	this->id = rhs.id;
}


template<class Model>
Cluster<Model>& Cluster<Model>::operator=(const Cluster<Model>& rhs) {
	if (this != &rhs){
		this->age = rhs.age;
		this->w = rhs.w;
		this->gamma = rhs.gamma;
		this->prm = rhs.prm;
		this->clusData = rhs.clusData;
		this->id = rhs.id;
	}
	return *this;
}


template<class Model>
void Cluster<Model>::updatePrm(){
	this->prm.update(this->clusData.begin(), this->clusData.end(), this->oldprm);
}

template<class Model>
void Cluster<Model>::finalize(){
	if(this->isEmpty()){
		this->age++;
	} else {
		this->oldprm.updateOld(this->clusData.begin(), this->clusData.end(), this->oldprm);
		this->age = 1;
	}
	this->clusData.clear();
}

template<class Model>
std::vector<uint64_t> Cluster<Model>::getAssignedIds() const{
	std::vector<uint64_t> asids;
	for (auto it = this->clusData.begin(); it != this->clusData.end(); ++it){
		asids.push_back(it->first);
	}
	return asids;
}

template<class Model>
void Cluster<Model>::assignData(uint64_t did, Model::Data& d){
	if (this->clusData.find(did) != this->clusData.end()){
		throw DataAlreadyInClusterException(did);
	}
	this->clusData[did] = d;
}

template<class Model>
Model::Data Cluster<Model>::deassignData(uint64_t did){
	if (this->clusData.find(did) == this->clusData.end()){
		throw DataNotInClusterException(did); 
	}
	Model::Data d = this->clusData[did];
	this->clusData.erase(did);
	return d;
}

template<class Model>
void Cluster<Model>::clearData(){
	this->clusData.clear();
}

template<class Model>
void Cluster<Model>::setID(uint64_t id){
	if (this->id >= 0){
		throw IDAlreadySetException(this->id, id);
	}
	this->id = id;
}

template<class Model>
double Cluster<Model>::distTo(const Model::Data& d) const{
	if (this->isEmpty()){
		throw ClusterEmptyDistanceException();
	}
	return this->prm.distTo(d);
}

template<class Model>
double Cluster<Model>::distToOld(const Model::Data& d) const{
	return this->prm.distToOld(d);
}

template<class Model>
double Cluster<Model>::cost(double lambda, double Q) const{
	return this->isEmpty() ? 0.0 :
		(this->age == 0 ? lambda : Q*this->age)
		+this->prm.cost(this->clusData.begin(), this->clusData.end(), this->gamma);
}

template<class Model>
bool Cluster<Model>::isEmpty() const{
	return this->clusData.empty();
}

template<class Model>
bool Cluster<Model>::isNew() const{
	return this->age == 0;
}

template<class Model>
bool Cluster<Model>::isPermanentlyDead() const{
	return Model::isPermanentlyDead(this->age);
}

template<class Model>
const Model::Parameter& Cluster<Model>::getPrm() const{
	return this->prm;
}

#define __CLUSTER_IMPL_HPP
#endif /* __CLUSTER_IMPL_HPP */
