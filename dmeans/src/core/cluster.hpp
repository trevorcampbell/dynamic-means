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
		Cluster(const Cluster<D, P>& rhs);
		Cluster<D, P>& operator=(const Cluster<D, P>& rhs);
		void setID(uint64_t id);
		uint64_t getID() const;
		void assignData(uint64_t id, const D& d);
		std::vector<uint64_t> getAssignedIds() const;
		D deassignData(uint64_t did);
		void clearData();
		void finalize();
		bool isEmpty() const;
		bool isNew() const;
		uint64_t getAge() const;
		P& getPrmRef();
		P& getOldPrmRef();
		const P& getPrm() const;
		const P& getOldPrm() const;
		typename std::map<uint64_t, D>::const_iterator data_cbegin() const;
		typename std::map<uint64_t, D>::const_iterator data_cend() const;
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
		class IDUnsetException{
			public:
				IDUnsetException(){
					std::cout << "Trying to get cluster with unset id." << std::endl;
				}
		};
};

#include "cluster_impl.hpp"
}
#define __CLUSTER_HPP
#endif /* __CLUSTER_HPP */
