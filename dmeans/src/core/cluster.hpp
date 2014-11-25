#ifndef __CLUSTER_HPP
#include<map>
#include<vector>
#include<iostream>
#include<iterator>
#include "../util/config.hpp"

namespace dmeans{
template <class D, class P>
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
		void finalize();
		double compareTo(const D& d) const;
		P& getPrm() const;
		P& getOldPrm() const;
	private:
		uint64_t id;
		uint64_t age;
		P prm, oldprm;
		std::map<uint64_t, D> clusData;

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
