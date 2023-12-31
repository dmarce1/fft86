/*
 * timer.hpp
 *
 *  Created on: Aug 23, 2023
 *      Author: dmarce1
 */

#ifndef TIMER_HPP_
#define TIMER_HPP_

#include <chrono>

class timer {
	std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
	double time;
public:
	inline timer() {
		time = 0.0;
	}
	inline void stop() {
		std::chrono::time_point<std::chrono::high_resolution_clock> stop_time = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> dur = stop_time - start_time;
		time += dur.count();
	}
	inline void start() {
		start_time = std::chrono::high_resolution_clock::now();
	}
	inline void reset() {
		time = 0.0;
	}
	inline double read() {
		return time;
	}
};


#endif /* TIMER_HPP_ */
