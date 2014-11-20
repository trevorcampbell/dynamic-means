#ifndef __EXPKERNEL_HPP
#include <Eigen/Dense>

class ExpKernelData{
	public:
		Eigen::VectorXd v;
		double distTo(const VectorData& rhs){
			return (this->v - rhs.v).squaredNorm();
		}
};

typedef std::map<int, Data<VectorData> >::iterator vector_dmap_iterator;
class ExpKernelParameter{
	public:
		Eigen::VectorXd v, vOld;
		void update(vector_dmap_iterator be, vector_dmap_iterator en, double gamma){
			v = gamma*vOld;
			double wt = gamma;
			for(auto it = be; it != en; ++it){
				v += it->second.d.v;
				wt += 1.0;
			}
			v /= wt;
		}
		double cost(vector_dmap_iterator be, vector_dmap_iterator en, double gamma){
			double c = gamma*(v-vOld).squaredNorm();
			for(auto it = be; it != en; ++it){
				c += (v-it->second.d.v).squaredNorm();
			}
			return c;
		}
		std::vector<uint64_t> updateOld(vector_dmap_iterator be, vector_dmap_iterator en, double gamma){
			Eigen::VectorXd tmpv = gamma*vOld;
			double wt = gamma;
			for(auto it = be; it != en; ++it){
				tmpv += it->second.d.v;
				wt += 1.0;
			}
			vOld = tmpv / wt;
			return std::vector<uint64_t>();
		}

		double distTo(const VectorData& vec, bool isActive){
			return (isActive ? (vec.v-this->v).squaredNorm() : (vec.v-this->vOld).squaredNorm());
		}
};


#define __EXPKERNEL_HPP
#endif /* __EXPKERNEL_HPP */
