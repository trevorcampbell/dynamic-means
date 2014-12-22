#ifndef __MATCHINGSPECTRAL_HPP
#include<Eigen/Dense>
#include "../util/eigensolver.hpp"
#include "../util/maxmatching.hpp"
#include <cassert>


namespace dmeans{
template<class Model, bool monoCheck>
class _MatchingSpectral{
	using Clus = Cluster<typename Model::Data, typename Model::Parameter>;
	public:
		_MatchingSpectral(const Config& cfg);
		double cluster(const std::vector<typename Model::Data>& obs, std::vector<Clus>& clus, const Model& model) const;
	private:
		bool verbose;
		EigenSolver::Type solverType;
		uint64_t nEigs;
		Config cfg;
		uint64_t nProjectionRestarts;

		MXd getKernelMatUpper(const std::vector<typename Model::Data>& obs, const std::vector<Clus>& clus, const Model& model) const;
		void findClosestConstrained(const MXd& ZV, MXd& X) const;
		void findClosestRelaxed(const MXd& Z, const MXd& X, MXd& V) const; 
		void orthonormalize(MXd& V) const; 
		vector<uint64_t> getLblsFromIndicatorMat(const MXd& X) const;
		double getOldNewMatching(const std::vector<uint64_t>& lbls, const std::vector<typename Model::Data>& obs, 
				const std::vector<Clus>& clus, const Model& model, std::map<uint64_t, uint64_t>& lblmap);
		class MonotonicityViolationException{
			public:
				MonotonicityViolationException(double prevobj, double obj, std::string funcname){
					std::cout << "Monotonicity violated! Prevobj = " << prevobj << " obj = " << obj << " after calling " << funcname << std::endl;
				}
		};
};

template<class Model>
using MatchingSpectralWithMonotonicityChecks = _MatchingSpectral<Model, true>; 

template<class Model>
using MatchingSpectral = _MatchingSpectral<Model, false>;

#include "matchingspectral_impl.hpp"
}
#define __MATCHINGSPECTRAL_HPP
#endif /* __MATCHINGSPECTRAL_HPP */
