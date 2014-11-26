#ifndef __SPARSEVECTORAPPROX_HPP

#include<Eigen/Dense>
#include<vector>

typedef Eigen::MatrixXd MXd;
typedef Eigen::VectorXd VXd;

namespace dmeans{
class SparseVectorApproximation{
	public:
		SparseVectorApproximation(uint64_t K, double eps);
		void fromKernelMatrix(MXd m, VXd acoeffs);
		void getApprox(std::vector<uint64_t>& vecs, std::vector<double>& acoeffs);
	private:
		uint64_t K;
		double eps;
		std::vector<uint64_t> bvecs;
		std::vector<double> coeffs;

		void projectOntoOrthogonalSubspace(VXd& acolsums, MXd& m, std::vector<uint64_t>& idcs, std::vector<uint64_t>& srtidcs, const uint64_t kMax);
		void removeRowCol(MXd& m, uint64_t k);
		void removeElem(VXd& v, uint64_t k);
};

#include "sparsevectorapprox_impl.hpp"
}

#define __SPARSEVECTORAPPROX_HPP
#endif /* __SPARSEAPPROX_HPP */
