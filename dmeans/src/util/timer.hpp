#ifndef __TIMER_HPP
#include<chrono>
namespace dmeans{
class Timer{
	public:
		typedef std::chrono::high_resolution_clock hrc;
		hrc::time_point t0;
		double now_s(){
			return (std::chrono::duration_cast< std::chrono::duration<double> >(hrc::now().time_since_epoch())).count();
		}
		void start(){
			t0 = hrc::now();
		}
		double elapsed_ms(){
			return (std::chrono::duration_cast< std::chrono::duration<double> >(hrc::now()-t0)).count();
		}
};
}
#define __TIMER_HPP
#endif /* __TIMER_HPP */
