#ifndef __CLUSTER_HPP
#include<map>
#include<vector>
#include<iostream>
#include<iterator>
namespace dmeans{
template <class D, class P>
class Cluster {
	public:
		const uint64_t id;

		Cluster(double lambda, double Q, double tau);
		void updatePrm();
		double cost() const;
		double gamma() const;
		void assignData(uint64_t id, D& d);
		std::vector<uint64_t> getAssignedIds() const;
		D deassignData(uint64_t did);
		void finalize();
		double distTo(const D& d) const;
		bool isEmpty() const;
		bool isNew() const;
		P getPrm() const;
	private:
		static uint64_t nextId;
		uint64_t age;
		double w, lambda, Q, tau;
		P prm;
		std::map<uint64_t, D> clusData;

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

template<class D, class P> 
uint64_t Cluster<D, P>::nextId = 0;

#include "cluster_impl.hpp"
}
#define __CLUSTER_HPP
#endif /* __CLUSTER_HPP */
