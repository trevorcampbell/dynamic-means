#ifndef __ITERATIVE_HPP
#include<vector>
#include<iostream>
#include<random>
#include "../core/cluster.hpp"
#include "../util/config.hpp"
#include "../util/random.hpp"
#include<string>

namespace dmeans{
template<class Model, bool monoCheck>
class _Iterative{
	using Clus = Cluster<typename Model::Data, typename Model::Parameter>;
	public:
		_Iterative(const Config& cfg);
		double cluster(const std::vector<typename Model::Data>& obs, std::vector<Clus>& clus, const Model& model) const;
	private:
		bool verbose;
		Config cfg;

		void initialLabelling(const std::vector< typename Model::Data>& obs, std::vector<Clus>& clus, const Model& model) const;
		bool labelUpdate(std::vector<Clus>& clus, const Model& model) const;
		void parameterUpdate(std::vector<Clus>& clus, const Model& model) const;
		double computeCost(const std::vector<Clus>& clus, const Model& model) const;
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
