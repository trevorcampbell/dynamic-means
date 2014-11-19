#ifndef __DYNMEANS_HPP
#include<vector>
#include<iostream>
#include<algorithm>
#include<boost/static_assert.hpp>
#include<boost/function.hpp>
#include<boost/bind.hpp>
#include<sys/time.h>
#include <ctime>

template <class Vec>
class DynMeans{
	public:
		DynMeans(double lambda, double Q, double tau, bool verbose = false);
		~DynMeans();
		
		//initialize a new step and cluster
		void cluster(std::vector<Vec>& newobservations, int nRestarts, std::vector<int>& finalLabels, std::vector<Vec>& finalParams, double& finalObj, double& tTaken);
		//reset DDP chain
		void reset();
	private:
		double lambda, Q, tau;
		bool verbose;
		std::vector<Vec> observations;
		//during each step, constants which are information about the past steps
		//once each step is complete, these get updated
		int nextLbl;
		std::vector<Vec> oldprms;
		std::vector<int> oldprmlbls;
		std::vector<double> weights;
		std::vector<int> ages;

		//tools to help with kmeans
		std::vector<Vec> getObsInCluster(int idx, std::vector<int> lbls); 
		void assignObservations(std::vector<int> assgnOrdering, std::vector<int>& lbls, std::vector<int>& cnts, std::vector<Vec>& prms);
		double setParameters(std::vector<int>& lbls, std::vector<int>& cnts, std::vector<Vec>& prms);
		std::vector<int> updateState(std::vector<int> lbls, std::vector<int> cnts, std::vector<Vec> prms);
};
#include "dynmeans_impl.hpp"
#define __DYNMEANS_HPP
#endif /* __DYNMEANS_HPP */
