#ifndef __DOTKERNELMODEL_HPP

#include <Eigen/Dense>
#include "../core/cluster.hpp"
#include "../util/sparsevectorapprox.hpp"
namespace dmeans{
template <int n>
class DotProductKernelModel{
	public:
		double lambda, Q, tau, spEps;
		uint64_t spK;
		
		class InvalidParameterException{
			public:
				InvalidParameterException(double lambda, double Q, double tau, uint64_t spK, double spEps){
					std::cout << "Error - Invalid parameters for DotProductKernelModel!" << std::endl;
					std::cout << "Lambda: " << lambda << std::endl;
					std::cout << "Q: " << Q << std::endl;
					std::cout << "tau: " << tau << std::endl;
					std::cout << "spK: " << spK << std::endl;
					std::cout << "spEps: " << spEps << std::endl;
				}
		};

		DotProductKernelModel(Config cfg){
			this->lambda = cfg.get("lambda", Config::REQUIRED, -1.0);
			this->Q = cfg.get("Q", Config::REQUIRED, -1.0);
			this->tau = cfg.get("tau", Config::REQUIRED, -1.0);
			this->spK = cfg.get("sparseApproximationSize", Config::REQUIRED, -1);
			this->spEps = cfg.get("sparseApproximationErrorThreshold", Config::OPTIONAL, 1.0e-6);
			if (lambda <= 0 || Q <= 0 || tau <= 0 || spK <= 0 || spEps <= 0){
				throw InvalidParameterException(lambda, Q, tau, spK, spEps);
			}
		}
	

		class Data{
			public:
				Data(){v.setZero();}
				Eigen::Matrix<double, n, 1> v;
		};
		class Parameter{
			public:
				Parameter(){vs.clear(); w = kp2p = 0;}
				std::vector< Eigen::Matrix<double, n, 1> > vs;
				std::vector< double > coeffs;
				double w, kp2p;
		};


		void update(const std::vector<Data>& obs){
			//nothing to update, dot product model doesn't depend on current data
		}

		bool isClusterDead(double age) const{
			return Q*age > lambda;
		}

		typedef typename std::map<uint64_t, Data>::const_iterator dmap_iterator; 

		bool exceedsNewClusterCost(const Data& d, double cost) const{
			return cost > lambda;
		}

		double getEigenvalueLowerThreshold() const{
			return lambda;
		}

		double getOldPenalty(const Cluster<Data, Parameter>& c) const {
			return Q*c.getAge();
		}
		double getNewPenalty() const {
			return lambda;
		}


		double oldWeight(const Cluster<Data, Parameter>& c) const {
			if (c.getAge() == 0){ return 0.0;}
			return 1.0/(1.0/c.getOldPrm().w +tau*c.getAge()); 
		}

		double kernelDD(const Data& d1, const Data& d2) const{
			return d1.v.dot(d2.v);
		}

		double kernelDOldP(const Data& d, const Cluster<Data, Parameter>& c) const{
			double kern = 0.0;
			for (uint64_t i = 0; i < c.getOldPrm().vs.size(); i++){
				kern += c.getOldPrm().coeffs[i]*d.v.dot(c.getOldPrm().vs[i]);
			}
			return kern;
		}

		double kernelOldPOldP(const Cluster<Data, Parameter>& c) const{
			return c.getOldPrm().kp2p;
		}

		double kernelPP(const Cluster<Data, Parameter>& c) const{
			return c.getPrm().kp2p;
		}

		double clusterCost(const Cluster<Data, Parameter>& c) const{
			if (c.isEmpty()){
				return 0.0;
			}
			double cost = 0.0;
			double age = c.getAge();
			double gamma = oldWeight(c);
			uint64_t N = c.getAssignedIds().size();
			if (c.isNew()){
				cost += lambda;
			} else {
				cost += Q*age;
			}
			cost += -(gamma+N)*kernelPP(c) +gamma*kernelOldPOldP(c);
			for(auto it = c.data_cbegin(); it != c.data_cend(); ++it){
				cost += kernelDD(it->second, it->second);
			}
			return cost;
		}

		void updatePrm(Cluster<Data, Parameter>& c) const{
			if (c.isEmpty()){
				c.getPrmRef() = c.getOldPrmRef();
				return;
			} else {
				double gamma = oldWeight(c);
				uint64_t N = c.getAssignedIds().size();

				std::vector< Eigen::Matrix<double, n, 1> > fullvs;
				std::vector< double > fullcoeffs;
				for (uint64_t i=0; i < c.getOldPrm().vs.size(); i++) {
					fullvs.push_back( c.getOldPrm().vs[i] );
					fullcoeffs.push_back( c.getOldPrm().coeffs[i]*gamma/(gamma+N) );
				}
				for (auto it = c.data_cbegin(); it != c.data_cend(); ++it){
					fullvs.push_back(it->second.v);
					fullcoeffs.push_back(1.0/(gamma+N));
				}

				SparseVectorApproximation spa(spK, spEps);
				MXd kmat = MXd::Zero(fullvs.size(), fullvs.size());
				VXd coeffvec = VXd::Zero(fullvs.size());
				for (uint64_t i = 0; i < fullvs.size(); i++){
					for (uint64_t j = 0; j <= i; j++){
						kmat(i, j) = kmat(j, i) = fullvs[i].dot(fullvs[j]);
					}
					coeffvec(i) = fullcoeffs[i];
				}
				spa.fromKernelMatrix(kmat, coeffvec);
				std::vector<uint64_t> approxvecs;
				std::vector<double> approxcoeffs;
				spa.getApprox(approxvecs, approxcoeffs);

				c.getPrmRef().vs.clear();
				c.getPrmRef().coeffs.clear();
				for(uint64_t i = 0; i < approxvecs.size(); i++){
					c.getPrmRef().vs.push_back(fullvs[approxvecs[i]]);
					c.getPrmRef().coeffs.push_back(approxcoeffs[i]);
				}

				c.getPrmRef().w = gamma+N;


				//now cache the kernel from prm->prm and prm->oldprm
				double kp2p = 0.0;
				for (uint64_t i = 0; i < c.getPrm().vs.size(); i++){
					for (uint64_t j = 0; j < c.getPrm().vs.size(); j++){
						kp2p += c.getPrm().coeffs[i]*c.getPrm().coeffs[j]*c.getPrm().vs[j].dot(c.getPrm().vs[i]);
					}
				}
				c.getPrmRef().kp2p = kp2p;
			}
		}

		void kernelAdd(Cluster<Data, Parameter>& c, const Data& d) const {
			double gamma = oldWeight(c);

			if (c.isEmpty()){
				c.getPrmRef().kp2p = kernelOldPOldP(c);
			}
			uint64_t N = c.getAssignedIds().size();
			double kp2p = c.getPrmRef().kp2p;
			kp2p *= (gamma+N)*(gamma+N);
			kp2p += kernelDD(d, d) + 2*gamma*kernelDOldP(d, c);
			for(auto it = c.data_cbegin(); it != c.data_cend(); ++it){
				kp2p += 2*kernelDD(d, it->second);
			}
			c.getPrmRef().kp2p = kp2p/((gamma+N+1)*(gamma+N+1));
		}

		void kernelSubtract(Cluster<Data, Parameter>& c, const Data& d) const {
			if (c.isEmpty()){
				c.getPrmRef().kp2p = 0.0;
				return;
			}
			double gamma = oldWeight(c);
			uint64_t N = c.getAssignedIds().size() + 1;
			double kp2p = c.getPrmRef().kp2p;
			kp2p *= (gamma+N)*(gamma+N);
			kp2p -= kernelDD(d, d) + 2*gamma*kernelDOldP(d, c);
			for(auto it = c.data_cbegin(); it != c.data_cend(); ++it){
				kp2p -= 2*kernelDD(d, it->second);
			}
			c.getPrmRef().kp2p = kp2p/((gamma+N-1)*(gamma+N-1));
		}

		double compare(const Cluster<Data, Parameter>& c, const Data& d) const{
			if (c.isEmpty()){
				double age = c.getAge();
				double gamma = oldWeight(c);
				return Q*age+gamma/(gamma+1.0)*(kernelDD(d, d) - 2*kernelDOldP(d, c) + kernelOldPOldP(c));
			} else {
				double gamma = oldWeight(c);
				uint64_t N = c.getAssignedIds().size();
				double ret = kernelDD(d, d) - 2*gamma*kernelDOldP(d, c)/(gamma+N);
				for(auto it = c.data_cbegin(); it != c.data_cend(); ++it){
					ret -= 2*kernelDD(d, it->second)/(gamma+N);
				}
				ret += kernelPP(c);
				return (gamma+N)/(gamma+N+1)*ret;
			}
		}

};
}


#define __DOTKERNELMODEL_HPP
#endif /* __DOTKERNELMODEL_HPP */
