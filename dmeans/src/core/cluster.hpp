#ifndef __CLUSTER_HPP
#include<map>
#include<vector>
#include<iostream>
#include<iterator>

namespace dmeans{
template <class Model>
class Cluster {
	public:
		Cluster();
		Cluster(const Cluster<Model>& rhs);
		Cluster<Model>& operator=(const Cluster<Model>& rhs);
		void setID(uint64_t id);
		void assignData(uint64_t id, typename Model::Data& d);
		std::vector<uint64_t> getAssignedIds() const;
		typename Model::Data deassignData(uint64_t did);
		void clearData();
		bool isEmpty() const;
		bool isNew() const;
		bool isPermanentlyDead() const;
		void finalize();
		void updatePrm();
		double cost() const;
		double compareTo(typename Model::Data& d) const;
		const typename Model::Parameter& getPrm() const;
	private:
		uint64_t id;
		uint64_t age;
		Model::Parameter prm, oldprm;
		std::map<uint64_t, typename Model::Data> clusData;

		class DataNotInClusterException{
			public:
				DataNotInClusterException(uint64_t did){
					std::cout << "Cluster does not own datapoint " << did << std::endl;
				}
		};
		class DataAlreadyInClusterException{
			public:
				DataAlreadyInClusterException(uint64_t did){
					std::cout << "Cluster already owns datapoint " << did << std::endl;
				}
		};

		class ClusterEmptyDistanceException{
			public:
				ClusterEmptyDistanceException(){
					std::cout << "Requesting distance to the new parameter of an empty cluster " << std::endl;
				}
		};
		class IDAlreadySetException{
			public:
				IDAlreadySetException(uint64_t cid, uint64_t nid){
					std::cout << "Trying to change cluster id " << cid << " to " << nid << std::endl;
				}
		};
};

#include "cluster_impl.hpp"
}
#define __CLUSTER_HPP
#endif /* __CLUSTER_HPP */
