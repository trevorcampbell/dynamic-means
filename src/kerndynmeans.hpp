#ifndef __KERNDYNMEANS_HPP
#include<vector>
#include<stack>
#include<queue>
#include<iostream>
#include<algorithm>
#include<limits>
#include<boost/static_assert.hpp>
#include<boost/function.hpp>
#include<boost/bind.hpp>
#include<sys/time.h>
#include <ctime>
#include "gurobi_c++.h" //note: the use of this library requires gurobi!
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

using namespace std;

typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SMXd;
typedef Eigen::Triplet<double> TD;
typedef Eigen::MatrixXd MXd;
typedef Eigen::VectorXd VXd;

template <class G>
class KernDynMeans{
	public:
		KernDynMeans(double lambda, double Q, double tau, bool verbose = false);
		~KernDynMeans();
		//initialize a new step and cluster
		void cluster(const G& aff, const int nRestarts, const int nCoarsest, std::vector<int>& finalLabels, double& finalObj, std::vector<double>& finalGammas, 
		std::vector<int>& finalPrmLbls, double& tTaken);
		//reset DDP chain
		void reset();
	private:
		//clusters a refinement level with kernelized dyn means batch updates
		//template so it works with C/D
		template <typename T> std::vector<int> clusterAtLevel(const T& aff, std::vector<int> lbls);
		//compute the dynamic means objective given the current labels
		template<typename T> double objective(const T& aff, const std::vector<int>& lbls) const;
		//compute a minimum weight bipartite matching
		std::map<int, int> getMinWtMatching(vector< pair<int, int> > nodePairs, vector<double> edgeWeights ) const;
		//get the minimum weight old/new cluster correspondence
		template <typename T> std::vector<int> updateOldNewCorrespondence(const T& aff, std::vector<int> lbls) const;
		//get the updated data labels via dyn means iteration
		template <typename T> std::vector<int> updateLabels(const T& aff, std::vector<int> lbls) const;
		//update the state after all iterations are done
		void finalizeStep(const G& aff, const vector<int>& lbls, vector<double>& prevgammas_out, vector<int>& prmlbls_out);
		//do a base clustering using spectral methods + minimum weight matching
		template <typename T> std::vector<int> baseCluster(const T& aff) const;
		//utility function to orthonormalize a square matrix
		void orthonormalize(MXd& V) const;
		void initializeSigma(const G& aff);

		GRBEnv* grbenv;
		double lambda, Q, tau;
		bool verbose;
		double sigma, sigmaUB, sigmaLB;//correction used to enforce positive definiteness

		//during each step, constants which are information about the past steps
		//once each step is complete, these get updated
		std::vector<int> oldprmlbls;// require a mapping from old parameters to old labels
		int maxLblPrevUsed;         //because clusters that die are removed entirely to save computation
		std::vector<double> weights;
		std::vector<int> ages;
		std::vector<double> agecosts;
		std::vector<double> gammas;
};

template <class G>
class CoarseGraph{ //only stores upper triangular information (multiplies by 2 when necessary)
	public:
		//function that constructs the coarse graph
		template <typename T> void coarsify(const T& aff);
		//similarity functions
		double diagSelfSimDD(const int i) const;
		double offDiagSelfSimDD(const int i) const;
		double selfSimPP(const int i) const;
		double simDD(const int i, const int j) const;
		double simDP(const int i, const int j) const;
		int getNodeCt(const int i) const;
		//input labels for this coarsified graph, get the labels for the original refined graph
		std::vector<int> getRefinedLabels(const std::vector<int>& lbls) const;
		//get the number of graph nodes
		int getNNodes() const;
		int getNOldPrms() const;
	private:
		int nOldPrms;
		std::vector< std::pair<int, int> > refineMap;
		std::vector<int> nodeCts;
		SMXd affdd, affdp;
		VXd daffdd, odaffdd, affpp;
};

////try to split a cluster
//template <typename T> std::vector<int> clusterSplit(std::vector<T>& data, std::vector<int> lbls);
////try to merge two clusters
//template <typename T> std::vector<int> clusterMerge(std::vector<T>& data, std::vector<int> lbls);

#include "kerndynmeans_impl.hpp"
#define __KERNDYNMEANS_HPP
#endif /* __DYNMEANS_HPP */
