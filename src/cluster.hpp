#ifndef __CLUSTER_HPP
#include<map>
#include<vector>
#include<iostream>
#include<iterator>
#include "data.hpp"

template <class P, class D>
class Cluster {
	public:
		const uint64_t id;

		Cluster();
		void updatePrm();
		double cost() const;
		double gamma() const;
		void assignData(Data& d);
		std::vector<uint64_t> getAssignedIds() const;
		Data<D> deassignData(uint64_t did);
		std::vector<uint64_t> finalize();
		double distTo(const Data<D>& d) const;
		bool isEmpty() const;
		bool isNew() const;
		P getPrm() const;
	private:
		static uint64_t nextId;
		uint64_t age;
		double w;
		P prm;
		std::map<uint64_t, Data<D> > clusData;

		class DataNotInClusterException{
			public:
				DataNotInClusterException(uint64_t cid, uint64_t did){
					std::cout << "Datapoint " << did << "is not in cluster " << cid << std::endl;
				}
		};
		class DataAlreadyInClusterException{
			public:
				DataAlreadyInClusterException(uint64_t cid, uint64_t did){
					std::cout << "Datapoint " << did << "is already in cluster " << cid << std::endl;
				}
		};
};
uint64_t template<class P> Cluster<P>::nextId = 0;

#include "cluster_impl.hpp"
#define __CLUSTER_HPP
#endif /* __CLUSTER_HPP */
