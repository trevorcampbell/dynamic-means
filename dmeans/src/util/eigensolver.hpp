#ifndef __EIGENSOLVER_HPP
#include<Eigen/Dense>
#include<iostream>
#include<random>

typedef Eigen::MatrixXd MXd;
typedef Eigen::VectorXd VXd;

template<class T>
class EigenSolver{
	public:
		enum Type{
			EIGEN_SELF_ADJOINT,
			REDSVD
		};

		EigenSolver(MXd& A_UpperTriangle, Type t, uint64_t nEigs = 0, double lowerThresh=-1);
		void redsvdSolver(SMXd& AUp, int r);
		void selfadjointSolver(MXd& A_Up);
		void getResults(VXd& eigvals, MXd& eigvecs);
	private:
		MXd eigvecs;
		VXd eigvals;
		void columnGramSchmidt(MXd& m);
		void pruneSmallEigs(double thresh);
		class MatrixNotSquareException{
			public:
				MatrixNotSquareException(uint64_t rows, uint64_t cols){
					std::cout << "Matrix input to EigenSolver nonsquare: rows = " << rows << " cols = " << cols << std::endl;
				}
		};
};

#include "eigensolver_impl.hpp"

#define __EIGENSOLVER_HPP
#endif /* __EIGENSOLVER_HPP */
