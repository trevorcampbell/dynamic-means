#ifndef __MSTKERNELMODEL_HPP
#include <Eigen/Dense>
#include <algorithm>
#include <fstream>
#include <string>
#include "../core/cluster.hpp"
#include "../util/sparsevectorapprox.hpp"
namespace dmeans{
template <int n>
class MSTKernelModel{
	public:
		double lambda, Q, tau, jth, spEps;
		uint64_t spK;
		class MST{
			public:
				std::vector< Eigen::Matrix<double, n, 1> > nodes;
				std::vector< std::vector<int> > parentLists;
				void construct(std::vector< Eigen::Matrix<double, n, 1> > data){
					nodes.clear(); parentLists.clear();
					nodes.push_back(data[0]);
					parentLists.push_back(std::vector<int>());

					std::vector<int> toAsgn(data.size()-1);
					std::iota(toAsgn.begin(), toAsgn.end(), 1);
					while(!toAsgn.empty()){
						double minCost = std::numeric_limits<double>::infinity();
						int minInd = -1;
						int minParent = -1;
						for(uint64_t i = 0; i < toAsgn.size(); i++){
							for (uint64_t j = 0; j < nodes.size(); j++){
								double di = (nodes[j] - data[toAsgn[i]]).norm();
								if (di < minCost){
									minCost = di;
									minInd = i;
									minParent = j;
								}
							}
						}
						nodes.push_back(data[toAsgn[minInd]]);
						parentLists.push_back(parentLists[minParent]);
						parentLists.back().push_back(minParent);
						toAsgn.erase(toAsgn.begin()+minInd);
					}
				}

				double thresh(double di, double jth) const{
					if (di < jth){ return 0.0; }
					else { return di; }
				}
				double dist(Eigen::Matrix<double, n, 1> u, Eigen::Matrix<double, n, 1> v, double jth) const{
					//get the closest node to u, v
					if ( (u-v).norm() < jth){return 0.0;}
					double minU, minV;
					int paU, paV;
					minU = minV = std::numeric_limits<double>::infinity();
					for(uint64_t i = 0; i < nodes.size(); i++){
						if ( (nodes[i]-u).norm() < minU){
							minU = (nodes[i]-u).norm();
							paU = i;
						}
						if ( (nodes[i]-v).norm() < minV){
							minV = (nodes[i]-v).norm();
							paV = i;
						}
					}
					//now get the path between U and V
					std::vector<int> pasU = parentLists[paU];
					std::vector<int> pasV = parentLists[paV];
					pasU.push_back(paU);
					pasV.push_back(paV);
					//find the first index where the two lists differ
					int firstDifferId = -1;
					for(uint64_t i =0 ; i < std::min(pasU.size(), pasV.size()); i++){
						if (pasU[i] != pasV[i]){
							firstDifferId = i;
							break;
						}
					}
					double di = thresh(minU, jth)*thresh(minU, jth)+thresh(minV, jth)*thresh(minV, jth);
					if (firstDifferId == -1){
						//the shorter list is a subset of the longer list
						if (pasU.size() < pasV.size()){
							for (uint64_t i = pasU.size(); i < pasV.size(); i++){
								double tmp = thresh((nodes[pasV[i]] - nodes[pasV[i-1]]).norm(), jth);
								di += tmp*tmp;
							}
						} else {
							for (uint64_t i = pasV.size(); i < pasU.size(); i++){
								double tmp = thresh((nodes[pasU[i]] - nodes[pasU[i-1]]).norm(), jth);
								di +=tmp*tmp;
							}
						}
					} else {
						//it isn't a subset
						for(uint64_t i = pasU.size()-1; i > (uint32_t)firstDifferId; i--){ //the cast is safe because we already know firstDifferId >= 0
							double tmp = thresh( (nodes[pasU[i-1]] - nodes[pasU[i]]).norm() , jth);
							di += tmp*tmp;
						}
						for(uint64_t i = pasV.size()-1; i > (uint32_t)firstDifferId; i--){
							double tmp = thresh( (nodes[pasV[i-1]] - nodes[pasV[i]]).norm() , jth);
							di += tmp*tmp;
						}
					}
					return sqrt(di);
				}

				void write(std::string fname, double jth){
					std::ofstream treeout(fname.c_str(), std::ios_base::trunc);
					treeout << nodes.size() << std::endl;
					for (uint64_t i = 0; i < nodes.size(); i++){
						if (parentLists[i].size() == 0){
							treeout << -1 << " " << dist(nodes[i], nodes[0], jth) << " " <<  nodes[i].transpose() << std::endl;
						} else {
							treeout << parentLists[i].back() << " " <<  dist(nodes[i], nodes[0], jth) << " " << nodes[i].transpose() << std::endl;
						}
					}
					treeout.close();
				}
		} datatree;

		MSTKernelModel(Config cfg){
			this->lambda = cfg.get("lambda", Config::REQUIRED, -1.0);
			this->Q = cfg.get("Q", Config::REQUIRED, -1.0);
			this->tau = cfg.get("tau", Config::REQUIRED, -1.0);
			this->jth = cfg.get("jumpThreshold", Config::REQUIRED, -1.0);
			this->spK = cfg.get("sparseApproximationSize", Config::REQUIRED, -1);
			this->spEps = cfg.get("sparseApproximationErrorThreshold", Config::OPTIONAL, 1.0e-6);
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
			//MST model does depend on the current data (it computes 
			//distances using the MST created by the current dataset)
			std::vector< Eigen::Matrix<double, n, 1> > obsd;
			for (uint64_t i = 0; i < obs.size(); i++){
				obsd.push_back(obs[i].v);
			}
			datatree.construct(obsd);
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
			double di = datatree.dist(d1.v, d2.v, jth);
			return exp(-di*di/(2*jth*jth));
		}

		double kernelDOldP(const Data& d, const Cluster<Data, Parameter>& c) const{
			double kern = 0.0;
			for (uint64_t i = 0; i < c.getOldPrm().vs.size(); i++){
				double di = datatree.dist(d.v, c.getOldPrm().vs[i], jth);
				kern += c.getOldPrm().coeffs[i]*exp(-di*di/(2*jth*jth));
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
						double di = datatree.dist(fullvs[i], fullvs[j], jth);
						kmat(i, j) = kmat(j, i) = exp(-di*di/(2*jth*jth));
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
						double di = datatree.dist(c.getPrm().vs[j], c.getPrm().vs[i], jth);
						kp2p += c.getPrm().coeffs[i]*c.getPrm().coeffs[j]*exp(-di*di/(2*jth*jth));
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


#define __MSTKERNELMODEL_HPP
#endif /* __MSTKERNELMODEL_HPP */
