#ifndef __SPARSEVECTORAPPROX_HPP

#include<Eigen/Dense>
#include<vector>

typedef Eigen::MatrixXd MXd;
typedef Eigen::VectorXd VXd;

class SparseVectorApproximation{
	public:
		SparseVectorApproximation(int K, double eps);
		void fromKernelMatrix(MXd m, VXd acoeffs);
		void getApprox(std::vector<int>& vecs, std::vector<double>& acoeffs);
	private:
		int K;
		double eps;
		std::vector<int> bvecs;
		std::vector<double> coeffs;

		void projectOntoOrthogonalSubspace(VXd& acolsums, MXd& m, std::vector<int>& idcs, std::vector<int>& srtidcs, const int kMax);
		void removeRowCol(MXd& m, int k);
		void removeElem(VXd& v, int k);
};

#include "sparsevectorapprox_impl.hpp"

#define __SPARSEVECTORAPPROX_HPP
#endif /* __SPARSEAPPROX_HPP */
