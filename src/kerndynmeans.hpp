#ifndef __KERNDYNMEANS_HPP
#include<vector>
#include<stack>
#include<iostream>
#include<algorithm>
#include<limits>
#include<boost/static_assert.hpp>
#include<boost/function.hpp>
#include<boost/bind.hpp>
#include<sys/time.h>
#include <ctime>
#include<dynmeans/specdynmeans.hpp>
#include "gurobi_c++.h" //note: the use of this library requires gurobi!

using namespace std;

typedef Eigen::MatrixXd MXd;
typedef Eigen::VectorXd VXd;
typedef Eigen::SparseMatrix<double> SMXd;
typedef Eigen::Triplet<double> TD;


//TODO
//TODO
//TODO
//Make numInClus reflect coarsification levels, i.e. datapoints should be weighted heavier if they
//subsume more nodes in the original graph
//TODO
//TODO
//TODO

template <class D, class P>
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
		double objective(std::vector<D>& data, std::vector<int> lbls);
		//compute a minimum weight bipartite matching
		std::map<int, int> getMinWtMatching(vector< pair<int, int> > nodePairs, vector<double> edgeWeights ) const;
		//get the minimum weight old/new cluster correspondence
		template <typename T> std::vector<int> updateOldNewCorrespondence(std::vector<T>& data, std::vector<int> lbls);
		//get the updated data labels
		template <typename T> std::vector<int> updateLabels(std::vector<T>& data, std::vector<int> lbls);
		//try to split a cluster
		template <typename T> std::vector<int> clusterSplit(std::vector<T>& data, std::vector<int> lbls);
		//try to merge two clusters
		template <typename T> std::vector<int> clusterMerge(std::vector<T>& data, std::vector<int> lbls);


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

//FROM SDM
		////during each step, constants which are information about the past steps
		//once each step is complete, these get updated
		////vector<double> weights;
		////vector<int> ages;
		////vector<double> gammas, agecosts;
		////vector<P> oldprms;
		////vector<int> oldprmlbls; // require a mapping from old parameters to old labels
		////						//because clusters that die are removed entirely to save computation
		////int nextlbl; //stores the next (unique) label to use for a new cluster

		//////spectral functions
		////void getKernelMat(const vector<D>& data, SMXd& kUpper);
		////void solveEigensystem(SMXd& kUpper, const int nEigs, EigenSolverType type, MXd& eigenvectors);
		////void findClosestConstrained(const MXd& ZV, MXd& X) const;
		////map<int, int> getOldNewMatching(vector< pair<int, int> > nodePairs, vector<double> edgeWeights ) const;
		////void findClosestRelaxed(const MXd& Z, const MXd& X, MXd& V) const; 
		////void orthonormalize(MXd& V) const; 
		////double getNormalizedCutsObj(const SMXd& spmatUpper, const vector<int>& lbls) const;
		//////helper functions
		////vector<int> getLblsFromIndicatorMat(const MXd& X) const;
		//////Update the state for the new batch of data (timestep the ddp)
		////void updateState(const vector<D>& data, const vector<int>& lbls);
		////std::tuple<MXd, VXd> redsvdEigenSolver(SMXd& AUp, int r);
		////void gramschmidt(MXd& m);

////spectral data class 
//class SD{
//	public:
//		V2d v;
//		//similarity function from data->data is just exp(-|| ||^2 / w^2)
//		double sim(const SD& rhs) const{
//			return exp(-(this->v-rhs.v).squaredNorm()/(2*0.1*0.1));
//		}
//};
//
////spectral parameter class 
//class SP{
//	public:
//		V2d v;
//		//constructs a new parameter at the mean of a set of vectors
//		SP(const vector<SD>& rhs){
//			this->v = V2d::Zero();
//			for (int i = 0; i < rhs.size(); i++){
//				this->v += rhs[i].v;
//			}
//			this->v /= rhs.size();
//		}
//		//similarity function from parameter->data is just exp(-|| ||^2 / w^2)
//		double sim(const SD& rhs) const{
//			return exp(-(this->v-rhs.v).squaredNorm()/(2*0.1*0.1));
//		}
//		//This is the typical prior-weighted least squares Dynamic Means paramter update
//		void update(const vector<SD>& rhs, const double gamma){
//			this->v *= gamma;
//			for (int i = 0; i < rhs.size(); i++){
//				this->v += rhs[i].v;
//			}
//			this->v /= (gamma+rhs.size());
//		}
//};





#include "kerndynmeans_impl.hpp"
#define __KERNDYNMEANS_HPP
#endif /* __DYNMEANS_HPP */
