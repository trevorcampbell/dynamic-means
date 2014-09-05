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

using namespace std;

typedef Eigen::MatrixXd MXd;
typedef Eigen::VectorXd VXd;
typedef Eigen::SparseMatrix<double> SMXd;
typedef Eigen::Triplet<double> TD;

template <class D, class C, class P>
class KernDynMeans{
	public:
		KernDynMeans(double lambda, double Q, double tau, bool verbose = false);
		~KernDynMeans();
		//initialize a new step and cluster
		void cluster(std::vector<D>& data, const int nRestarts, const int nCoarsest, std::vector<int>& finalLabels, double& finalObj, double& tTaken);
		//reset DDP chain
		void reset();
	private:
		//clusters a refinement level with kernelized dyn means batch updates
		//tempalte so it works with C/D
		template <typename T> std::vector<int> clusterAtLevel(std::vector<T>& data, std::vector<int> initlbls);
		//expands the labels to the new refinement level
		std::vector<int> refine(std::vector< std::pair<int, int> > merges, std::vector<int> lbls);
		//template function for coarsify so it works with both C and D types
		template<typename T> std::pair< std::vector<C>, std::vector<std::pair<int, int> > > coarsify(std::vector<T>& data);
		//compute the dynamic means objective given the current labels
		template<typename T> double objective(std::vector<T>& data, std::vector<int> lbls);
		//compute a minimum weight bipartite matching
		std::map<int, int> getMinWtMatching(vector< pair<int, int> > nodePairs, vector<double> edgeWeights ) const;
		//get the minimum weight old/new cluster correspondence
		template <typename T> std::vector<int> updateOldNewCorrespondence(std::vector<T>& data, std::vector<int> lbls);
		//get the updated data labels
		template <typename T> std::vector<int> updateLabels(std::vector<T>& data, std::vector<int> lbls);
		////try to split a cluster
		//template <typename T> std::vector<int> clusterSplit(std::vector<T>& data, std::vector<int> lbls);
		////try to merge two clusters
		//template <typename T> std::vector<int> clusterMerge(std::vector<T>& data, std::vector<int> lbls);
		void updateState(const vector<D>& data, const vector<int>& lbls);
		std::vector<int> baseCluster(std::vector<T>& data);
		void orthonormalize(MXd& V) const; 


		GRBEnv* grbenv;
		double lambda, Q, tau;
		bool verbose;

		//during each step, constants which are information about the past steps
		//once each step is complete, these get updated
		std::vector<P> oldprms;
		std::vector<int> oldprmlbls;// require a mapping from old parameters to old labels
		int nextlbl;                //because clusters that die are removed entirely to save computation
		std::vector<double> weights;
		std::vector<int> ages;
		std::vector<double> agecosts;
		std::vector<double> gammas;
};


#include "kerndynmeans_impl.hpp"
#define __KERNDYNMEANS_HPP
#endif /* __DYNMEANS_HPP */
