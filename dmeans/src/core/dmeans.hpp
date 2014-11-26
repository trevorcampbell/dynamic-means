#ifndef __DMEANS_HPP
#include<vector>
#include<iostream>
#include "../util/timer.hpp"
#include "../util/results.hpp"
#include "cluster.hpp"
#include "../util/config.hpp"
#include "../util/random.hpp"

namespace dmeans{
template<class Model, template<typename> class Alg>
class DMeans{
	public:
		DMeans(Config cfg);
		//initialize a new step and cluster
		Results<Model> cluster(const std::vector<typename Model::Data>& obs);
		//reset DDP chain
		void reset();
	private:
		bool verbose;
		uint64_t nextLabel;
		uint64_t nRestarts;
		Config cfg;
		Model model;

		std::vector<Cluster<typename Model::Data, typename Model::Parameter> > clusters;
		Timer timer;

		Results<Model> getResults() const;
		void finalize();
		void restart();
		void labelNewClusters();
};

#include "dmeans_impl.hpp"

}
#define __DMEANS_HPP
#endif 
