#ifndef __ITERATIVE_DYNMEANS_HPP
#include<vector>
#include<iostream>
#include<random>
#include "../util/timer.hpp"
#include "../util/results.hpp"
#include "../model/data.hpp"
#include "../model/cluster.hpp"

template<class D, class P>
class IterativeDMeans{
	public:
		IterativeDMeans(double lambda, double Q, double tau, bool verbose = false, int seed = -1);
		//initialize a new step and cluster
		Results cluster(std::vector< Data<D> >& obs, uint64_t nRestarts);
		//reset DDP chain
		void reset();
	private:
		double lambda, Q, tau;
		bool verbose;

		std::vector< Cluster<D, P> > clusters;
		Timer timer;

		Results computeResults();
		void initialLabelling(std::vector< Data<D> >& obs);
		bool labelUpdate();
		void parameterUpdate();
		class MonotonicityViolationException{
			public:
				MonotonicityViolationException(double prevobj, double obj, string funcname){
					std::cout << "Monotonicity violated! Prevobj = " << prevobj << " obj = " << obj << " after calling " << funcname << std::endl;
				}
		};
}

#include "iterative_dmeans_impl.hpp"

#define __ITERATIVE_DYNMEANS_HPP
#endif /* __ITERATIVE_DYNMEANS_HPP */
