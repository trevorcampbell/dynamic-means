#ifndef __DATA_HPP
#include <set>
#include <iostream>

template <class D>
class Data{
	public:
		const uint64_t id;
		const D d;

		Data(uint64_t _id, D& _d): id(_id), d(_d){
			if (usedIds.find(id) != usedIds.end()){
				throw IDAlreadyUsedException(id);
			} else {
				usedIds.insert(id);
			}
		}

		Data(const Data<D>& rhs) : id(rhs.id), d(rhs.d){
		}
		double distTo(const Data<D>& rhs) const{
			return this->d.distTo(rhs.d);
		}
	private:
		static std::set<uint64_t> usedIds;
		class IDAlreadyUsedException{
			public:
				IDAlreadyUsedException(uint64_t uid){
					std::cout << "ID " << uid << " already in use, cannot create a new datapoint with that ID" << std::endl;
				}
		};
};
template<class D> std::set<uint64_t> Data<D>::usedIds;
#define __DATA_HPP
#endif /* __DATA_HPP */
