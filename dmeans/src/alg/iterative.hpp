#ifndef __ITERATIVE_HPP
#include<vector>
#include<iostream>
#include<random>
#include "../core/cluster.hpp"

namespace dmeans{
template<class D, class P, bool M>
class _Iterative{
	public:
		void cluster(std::map<uint64_t, D>& obs, std::map<uint64_t, Cluster<D, P> >& clus, double lambda, double Q, double tau, bool verbose);
	private:
		double computeCost();
		void initialLabelling(std::map<uint64_t, D>& obs);
		bool labelUpdate();
		void parameterUpdate();
		class MonotonicityViolationException{
			public:
				MonotonicityViolationException(double prevobj, double obj, const char* funcname){
					std::cout << "Monotonicity violated! Prevobj = " << prevobj << " obj = " << obj << " after calling " << funcname << std::endl;
				}
		};
};

template<class D, class P>
using _Iterative<D, P, true>  = IterativeWithMonotonicityChecks<D, P>;

template<class D, class P>
using _Iterative<D, P, false> = Iterative<D, P>;

#include "iterative_dmeans_impl.hpp"

}
#define __ITERATIVE_HPP
#endif 
