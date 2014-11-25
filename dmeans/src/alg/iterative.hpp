#ifndef __ITERATIVE_HPP
#include<vector>
#include<iostream>
#include<random>
#include "../core/cluster.hpp"
#include "../util/config.hpp"
#include<string>

namespace dmeans{
template<class Model, bool monoCheck>
class _Iterative{
	public:
		_Iterative(Config cfg);
		double cluster(std::map<uint64_t, typename Model::Data>& obs, std::vector<Cluster<Model> >& clus);
	private:
		bool verbose;
		Config cfg;

		void initialLabelling(std::map<uint64_t, typename Model::Data>& obs, std::vector< Cluster<Model> >& clus);
		bool labelUpdate(std::vector< Cluster<Model> >& clus);
		void parameterUpdate(std::vector< Cluster<Model> >& clus);
		double computeCost(std::vector< Cluster<Model> >& clus);
		class MonotonicityViolationException{
			public:
				MonotonicityViolationException(double prevobj, double obj, std::string funcname){
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
