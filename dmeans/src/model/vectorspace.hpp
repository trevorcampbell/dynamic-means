#ifndef __VECTORSPACE_HPP
#include <Eigen/Dense>
namespace dmeans{
template <int n>
class VectorData{
	public:
		Eigen::Matrix<double, n, 1> v;
		double distTo(const VectorData<n>& rhs) const{
			return (this->v - rhs.v).squaredNorm();
		}
};

template <int n>
using vector_dmap_iterator = typename std::map<uint64_t, VectorData<n> >::const_iterator;

template <int n>
class VectorParameter{
	public:
		Eigen::Matrix<double, n, 1> v, vOld;
		VectorParameter(){
			v = vOld = Eigen::VectorXd::Zero(n);
		}
		void update(vector_dmap_iterator<n> be, vector_dmap_iterator<n> en, double gamma){
			v = gamma*vOld;
			double wt = gamma;
			for(auto it = be; it != en; ++it){
				v += it->second.v;
				wt += 1.0;
			}
			v /= wt;
		}
		void updateOld(vector_dmap_iterator<n> be, vector_dmap_iterator<n> en, double gamma){
			Eigen::VectorXd tmpv = gamma*vOld;
			double wt = gamma;
			for(auto it = be; it != en; ++it){
				tmpv += it->second.v;
				wt += 1.0;
			}
			vOld = tmpv / wt;
		}

		double distTo(const VectorData<n> & vec) const{
			return (vec.v-this->v).squaredNorm();
		}

		double distToOld(const VectorData<n> & vec) const{
			return (vec.v-this->vOld).squaredNorm();
		}

		double cost(vector_dmap_iterator<n> be, vector_dmap_iterator<n> en, double gamma) const{
			double cost = 0;
			cost += gamma*(v-vOld).squaredNorm();
			for(auto it = be; it != en; ++it){
				cost += (it->second.v-v).squaredNorm();
			}
			return cost;
		}
};
}
#define __VECTORSPACE_HPP
#endif /* __VECTORSPACE_HPP */
