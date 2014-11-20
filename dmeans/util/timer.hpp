#ifndef __TIMER_HPP
#include<chrono>

class Timer{
	typedef std::chrono::high_resolution_clock hrc;
	typedef std::chrono::milliseconds ms;
	hrc::time_point t0;
	double now_ms(){
		return std::chrono::duration_cast<ms>(hrc::now().time_since_epoch()).count();
	}
	void start(){
		t0 = hrc::now();
	}
	double elapsed_ms(){
		return std::chrono::duration_cast<ms>(hrc::now()-t0).count();
	}
};
#define __TIMER_HPP
#endif /* __TIMER_HPP */
