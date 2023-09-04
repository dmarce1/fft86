/*
 * util.hpp
 *
 *  Created on: Aug 23, 2023
 *      Author: dmarce1
 */

#ifndef UTIL123_HPP_
#define UTIL123_HPP_

#include <fft86/vec.hpp>

#include <stdlib.h>

inline constexpr int round_down(int x, int y) {
	return y * (x / y);
}

inline constexpr int round_up(int x, int y) {
	return y * ((x - 1) / y + 1);
}

inline double rand1() {
	return (rand() + 0.5) / (RAND_MAX + 1.0);
}


#endif /* UTIL_HPP_ */
