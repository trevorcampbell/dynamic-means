#ifndef __SPECTRAL_HPP
#include<Eigen/Dense>
#include "../util/eigensolver.hpp"


namespace dmeans{
template<class Model, bool monoCheck>
class Spectral{
	using Clus = Cluster<typename Model::Data, typename Model::Parameter>;
	public:
		Spectral(const Config& cfg);
		double cluster(const std::vector<typename Model::Data>& obs, std::vector<Clus>& clus, const Model& model) const;
	private:
		bool verbose;
		EigenSolver::Type solverType;
		uint64_t nEigs;
		Config cfg;
		uint64_t nProjectionRestarts;

		MXd& getKernelMatUpper(const std::vector<typename Model::Data>& obs, const std::vector<Clus>& clus, const Model& model) const;
		void findClosestConstrained(const MXd& ZV, MXd& X) const;
		void findClosestRelaxed(const MXd& Z, const MXd& X, MXd& V) const; 
		void orthonormalize(MXd& V) const; 
		double getNormalizedCutsObj(const MXd& mUp, const vector<int>& lbls) const;
		vector<int> getLblsFromIndicatorMat(const MXd& X) const;
};

#include "spectral_impl.hpp"
#define __SPECTRAL_HPP
#endif /* __SPECTRAL_HPP */
