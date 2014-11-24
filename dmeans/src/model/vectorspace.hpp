#ifndef __VECTORSPACE_HPP
#include <Eigen/Dense>
namespace dmeans{
class VectorData{
	public:
		Eigen::VectorXd v;
		double distTo(const VectorData& rhs) const{
			return (this->v - rhs.v).squaredNorm();
		}
};

typedef std::map<uint64_t, VectorData >::const_iterator vector_dmap_iterator;
class VectorParameter{
	public:
		Eigen::VectorXd v, vOld;
		void update(vector_dmap_iterator be, vector_dmap_iterator en, double gamma){
			v = gamma*vOld;
			double wt = gamma;
			for(auto it = be; it != en; ++it){
				v += it->second.v;
				wt += 1.0;
			}
			v /= wt;
		}
		double cost(vector_dmap_iterator be, vector_dmap_iterator en, double gamma) const{
			double c = gamma*(v-vOld).squaredNorm();
			for(auto it = be; it != en; ++it){
				c += (v-it->second.v).squaredNorm();
			}
			return c;
		}
		void updateOld(vector_dmap_iterator be, vector_dmap_iterator en, double gamma){
			Eigen::VectorXd tmpv = gamma*vOld;
			double wt = gamma;
			for(auto it = be; it != en; ++it){
				tmpv += it->second.v;
				wt += 1.0;
			}
			vOld = tmpv / wt;
		}

		double distTo(const VectorData & vec, bool isActive) const{
			return (isActive ? (vec.v-this->v).squaredNorm() : (vec.v-this->vOld).squaredNorm());
		}
};
}
#define __VECTORSPACE_HPP
#endif /* __VECTORSPACE_HPP */
