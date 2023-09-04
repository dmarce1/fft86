/*
 * twiddles.hpp
 *
 *  Created on: Aug 22, 2023
 *      Author: dmarce1
 */

#ifndef TWIDDLES_HPP_
#define TWIDDLES_HPP_

#include <complex>
#include <vector>
#include <memory>

const std::complex<double>* get_twiddles(int N);
const double* get_cos_twiddles(int N);
const double* get_sin_twiddles(int N);
const double* get_rotated_cos_twiddles(int R, int N);
const double* get_rotated_sin_twiddles(int R, int N);
std::shared_ptr<std::vector<std::vector<double>>> get_6step_cos_twiddles(int N1, int N2);
std::shared_ptr<std::vector<std::vector<double>>> get_6step_sin_twiddles(int N1, int N2);

template<int N>
class cos_twiddle_array {
	static constexpr int M = N;
	double W[M];
public:
	constexpr double operator[](int i) const {
		return W[i];
	}
	constexpr cos_twiddle_array() :
			W() {
		for (int n = 0; n < M; n++) {
			W[n] = std::cos(2.0 * M_PI * n / N);
		}
		W[0] = 1.0;
		if (N % 2 == 0) {
			W[N / 2] = -1.0;
		}
		if (N % 4 == 0) {
			W[N / 4] = 0.0;
			W[3 * N / 4] = 0.0;
		}
	}
};

template<int N>
class sin_twiddle_array {
	static constexpr int M = N;
	double W[M];
public:
	constexpr double operator[](int i) const {
		return W[i];
	}
	constexpr sin_twiddle_array() :
			W() {
		for (int n = 0; n < M; n++) {
			W[n] = std::sin(2.0 * M_PI * n / N);
		}
		W[0] = 0.0;
		if (N % 2 == 0) {
			W[N / 2] = 0.0;
		}
		if (N % 4 == 0) {
			W[N / 4] = 1.0;
			W[3 * N / 4] = -1.0;
		}
	}
};

#endif /* TWIDDLES_HPP_ */
