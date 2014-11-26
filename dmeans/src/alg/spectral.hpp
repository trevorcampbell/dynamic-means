#ifndef __SPECTRAL_HPP

namespace dmeans{
template<class Model, bool monoCheck>
class _Iterative{
	using Clus = Cluster<typename Model::Data, typename Model::Parameter>;
	public:
		_Iterative(const Config& cfg);
		double cluster(const std::vector<typename Model::Data>& obs, std::vector<Clus>& clus, const Model& model) const;
	private:
		bool verbose;
		Config cfg;

		void initialLabelling(const std::vector< typename Model::Data>& obs, std::vector<Clus>& clus, const Model& model) const;
		bool labelUpdate(std::vector<Clus>& clus, const Model& model) const;
		void parameterUpdate(std::vector<Clus>& clus, const Model& model) const;
		double computeCost(const std::vector<Clus>& clus, const Model& model) const;
		class MonotonicityViolationException{
			public:
				MonotonicityViolationException(double prevobj, double obj, std::string funcname){
					std::cout << "Monotonicity violated! Prevobj = " << prevobj << " obj = " << obj << " after calling " << funcname << std::endl;
				}
		};
};

#include "spectral_impl.hpp"
#define __SPECTRAL_HPP
#endif /* __SPECTRAL_HPP */
