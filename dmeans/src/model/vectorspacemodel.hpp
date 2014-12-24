#ifndef __VECTORSPACEMODEL_HPP
#include <Eigen/Dense>
#include "../core/cluster.hpp"
namespace dmeans{
template <int n>
class VectorSpaceModel{
	public:
		double lambda, Q, tau;
		class InvalidParameterException{
			public:
				InvalidParameterException(double lambda, double Q, double tau){
					std::cout << "Error - Invalid parameters for DotProductKernelModel!" << std::endl;
					std::cout << "Lambda: " << lambda << std::endl;
					std::cout << "Q: " << Q << std::endl;
					std::cout << "tau: " << tau << std::endl;
				}
		};


		VectorSpaceModel(Config cfg){
			this->lambda = cfg.get("lambda", Config::REQUIRED, -1.0);
			this->Q = cfg.get("Q", Config::REQUIRED, -1.0);
			this->tau = cfg.get("tau", Config::REQUIRED, -1.0);
			if (lambda <= 0 || Q <= 0 || tau <= 0){
				throw InvalidParameterException(lambda, Q, tau);
			}
		}
		class Data{
			public:
				Data(){v.setZero();}
				Eigen::Matrix<double, n, 1> v;
		};
		class Parameter{
			public:
				Parameter(){v.setZero(); w = 0;}
				Eigen::Matrix<double, n, 1> v;
				double w;
		};

		typedef typename std::map<uint64_t, Data>::const_iterator dmap_iterator; 

		void update(const std::vector<Data>& obs) {
			//nothing to do here, vector space model doesn't depend on current dataset
		}

		bool isClusterDead(double age) const{
			return Q*age > lambda;
		}

		bool exceedsNewClusterCost(const Data& d, double cost) const{
			return cost > lambda;
		}

		double clusterCost(const Cluster<Data, Parameter>& c) const{
			if (c.isEmpty()){
				return 0.0;
			}
			double cost = 0.0;
			if (c.isNew()){
				cost += lambda;
			} else {
				double age = c.getAge();
				double gamma = 1.0/(1.0/c.getOldPrm().w +tau*age); //cluster old penalty
				cost += Q*age+gamma*(c.getOldPrm().v-c.getPrm().v).squaredNorm();
			}
			for(auto it = c.data_cbegin(); it != c.data_cend(); ++it){
				cost += (it->second.v-c.getPrm().v).squaredNorm();
			}
			return cost;
		}

		void updatePrm(Cluster<Data, Parameter>& c) const{
			if (c.isEmpty()){
				c.getPrmRef() = c.getOldPrmRef();
			} else {
				double age = c.getAge();
				double gamma = 1.0/(1.0/c.getOldPrm().w +tau*age); //cluster old penalty
				c.getPrmRef().v = gamma*c.getOldPrm().v;
				double wt = gamma;
				for(auto it = c.data_cbegin(); it != c.data_cend(); ++it){
					c.getPrmRef().v += it->second.v;
					wt += 1.0;
				}
				c.getPrmRef().v /= wt;
				c.getPrmRef().w = wt;
			}
		}

		double compare(const Cluster<Data, Parameter>& c, const Data& d) const{
			double age = c.getAge();
			double gamma = 1.0/(1.0/c.getOldPrm().w +tau*age);
			return (!c.isEmpty()) ? (d.v - c.getPrm().v).squaredNorm() : Q*age+gamma/(gamma+1.0)*(d.v-c.getOldPrm().v).squaredNorm();
		}
};
}
#define __VECTORSPACEMODEL_HPP
#endif /* __VECTORSPACE_HPP */
