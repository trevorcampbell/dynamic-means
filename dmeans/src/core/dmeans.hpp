#ifndef __DMEANS_HPP
#include<vector>
#include<iostream>
#include<random>
#include "../util/timer.hpp"
#include "../util/results.hpp"
#include "cluster.hpp"

namespace dmeans{
template<class D, class P, template<typename, typename> class A>
class DMeans{
	public:
		DMeans(double lambda, double Q, double tau, bool verbose = false, int seed = -1);
		//initialize a new step and cluster
		Results<P> cluster(std::map<uint64_t, D>& obs, uint64_t nRestarts);
		//reset DDP chain
		void reset();
	private:
		double lambda, Q, tau;
		bool verbose;

		std::map< uint64_t, Cluster<D, P> > clusters;
		Timer timer;

		Results<P> computeResults() const;
		void finalize();
		void restart();
};

#include "dmeans_impl.hpp"

}
#define __DMEANS_HPP
#endif 
