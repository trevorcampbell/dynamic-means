#ifndef __KERNELIZED_HPP
#include<vector>
#include<iostream>
#include<random>
#include "../core/cluster.hpp"
#include "../util/config.hpp"
#include "../util/random.hpp"
#include<string>

namespace dmeans{
template<class Model, bool monoCheck>
class _Kernelized{
	using Clus = Cluster<typename Model::Data, typename Model::Parameter>;
	public:
		_Kernelized(const Config& cfg);
		double cluster(const std::vector<typename Model::Data>& obs, std::vector<Clus>& clus, const Model& model) const;
	private:
		bool verbose;

		void initialLabelling(const std::vector< typename Model::Data>& obs, std::vector<Clus>& clus, const Model& model) const;
		bool labelUpdate(std::vector<Clus>& clus, const Model& model) const;
		void parameterUpdate(std::vector<Clus>& clus, const Model& model) const;
		double computeCost(const std::vector<Clus>& clus, const Model& model) const;
		class MonotonicityViolationException{
			public:
				MonotonicityViolationException(double prevobj, double obj, std::string funcname){
					std::cout << "Monotonicity violated! Prevobj = " << prevobj << " obj = " << obj << " after calling " << funcname << std::endl;
					std::cout << "Diff: " << obj-prevobj << std::endl;
				}
		};
};

template<class Model>
using KernelizedWithMonotonicityChecks = _Kernelized<Model, true>; 

template<class Model>
using Kernelized = _Kernelized<Model, false>;

#include "kernelized_impl.hpp"

}
#define __KERNELIZED_HPP
#endif 
