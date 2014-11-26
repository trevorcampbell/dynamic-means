#ifndef __RANDOM_HPP
#include<random>
#include "timer.hpp"

//singleton design pattern for random stuff
//singletons ensure that all randomness in dynamic means
//comes from a single source 
//helps with repeatability during experiments
namespace dmeans{
class RNG{
	public:
		static std::mt19937& get(){
			static RNG instance;
			return instance.rng_;
		}
	private:
		std::mt19937 rng_;
		RNG(){
			Timer ti;
			rng_.seed(ti.now_ms())
		}
		RNG(const RNG&);
		void operator=(const RNG&);
};

}
#define __RANDOM_HPP
#endif /* __RANDOM_HPP */
