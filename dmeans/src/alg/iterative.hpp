#ifndef __ITERATIVE_HPP
#include<vector>
#include<iostream>
#include<random>
#include "../core/cluster.hpp"

namespace dmeans{
template<class Model, bool monoCheck>
class _Iterative{
	public:
		_Iterative(bool verbose);
		double cluster(std::map<uint64_t, typename Model::Data>& obs, std::map<uint64_t, Cluster<Model> >& clus, double lambda, double Q, double tau, bool verbose);
	private:
		bool verbose;

		void initialLabelling(std::map<uint64_t, typename Model::Data>& obs, std::map<uint64_t, Cluster<Model> >& clus, double lambda);
		bool labelUpdate(std::map<uint64_t, Cluster<Model> >& clus, double lambda);
		void parameterUpdate(std::map<uint64_t, Cluster<Model> >& clus);
		double computeCost(std::map<uint64_t, Cluster<Model> >& clus, double lambda, double Q);
		class MonotonicityViolationException{
			public:
				MonotonicityViolationException(double prevobj, double obj, const char* funcname){
					std::cout << "Monotonicity violated! Prevobj = " << prevobj << " obj = " << obj << " after calling " << funcname << std::endl;
				}
		};
};

template<class Model>
using IterativeWithMonotonicityChecks = _Iterative<Model, true>; 

template<class Model>
using Iterative = _Iterative<Model, false>;

#include "iterative_impl.hpp"

}
#define __ITERATIVE_HPP
#endif 
