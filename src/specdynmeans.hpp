#ifndef __SPECDYNMEANS_HPP
#include<map>
#include<vector>
#include<string>
#include<iostream>
#include<algorithm>
#include<sys/time.h>
#include <ctime>
#include <random>
#include "gurobi_c++.h" //note: the use of this library requires gurobi!
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

using namespace std;

typedef Eigen::MatrixXd MXd;
typedef Eigen::VectorXd VXd;
typedef Eigen::SparseMatrix<double> SMXd;
typedef Eigen::Triplet<double> TD;

//Spectral Dynamic Means
//Constructor accepts 3 algorithmic parameters: labmda, Q, and tau -- See T Campbell et al, "Dynamic Clustering via Asymptotics of the Dependent Dirichlet Process", NIPS 2013
//Intuitive description (all three parameters >= 0)
//Lambda (cluster birth): Larger -> fewer clusters
//Q (cluster death): Larger -> clusters die faster
//tau (cluster motion): Larger -> clusters assumed to be more mobile, i.e. capable of moving farther in each timestep

//Typical usage of this library:
//1) Construct the SpecDynMeans object
//2) Chop your data up into suitable time windows
//3) Run SpecDynMeans::cluster on each time window in sequence


template <typename G>
class SpecDynMeans{
	public:
		enum EigenSolverType{
			EIGEN_SELF_ADJOINT,
			REDSVD
		};
		SpecDynMeans(double lamb, double Q, double tau, bool verbose = false, int seed = -1);
		~SpecDynMeans();
		void cluster(const G& aff, const int nRestarts, const int nClusMax, EigenSolverType type, vector<int>& finalLabels, double& finalObj, 
				std::vector<double>& finalGammas, std::vector<int>& finalPrmLbls, double& tTaken);

		//reset DDP chain
		void reset();

	private:
		GRBEnv* grbenv;
		mt19937 rng;
		double lamb, Q, tau;
		bool verbose;

		//during each step, constants which are information about the past steps
		//once each step is complete, these get updated
		vector<double> weights;
		vector<int> ages;
		vector<double> gammas, agecosts;
		vector<int> oldprmlbls; // require a mapping from old parameters to old labels
								//because clusters that die are removed entirely to save computation
		int maxLblPrevUsed; //stores the next (unique) label to use for a new cluster

		//spectral functions
		void getKernelMat(const G& aff, SMXd& kUpper);
		void solveEigensystem(SMXd& kUpper, const int nEigs, EigenSolverType type, MXd& eigenvectors);
		void findClosestConstrained(const MXd& ZV, MXd& X) const;
		map<int, int> getOldNewMatching(vector< pair<int, int> > nodePairs, vector<double> edgeWeights ) const;
		void findClosestRelaxed(const MXd& Z, const MXd& X, MXd& V) const; 
		void orthonormalize(MXd& V) const; 
		double getNormalizedCutsObj(const SMXd& spmatUpper, const vector<int>& lbls) const;
		//helper functions
		vector<int> getLblsFromIndicatorMat(const MXd& X) const;
		//Update the state for the new batch of data (timestep the ddp)
		void finalizeStep(const G& aff, const vector<int>& lbls, vector<double>& prevgammas_out, vector<int>& prmlbls_out);
		std::tuple<MXd, VXd> redsvdEigenSolver(SMXd& AUp, int r);
		void gramschmidt(MXd& m);
};


#include "specdynmeans_impl.hpp"
#define __SPECDYNMEANS_HPP
#endif /* __SPECDYNMEANS_HPP */


