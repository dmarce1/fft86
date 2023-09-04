//============================================================================
// Name        : fft86.cpp
// Author      : Dominic Marcello
// Version     :
// Copyright   : Copyright 2023
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <fft86/fftw.hpp>
#include <fft86/fft.hpp>
#include <fft86/util.hpp>
#include <fft86/defs.hpp>

#include <algorithm>
#include <complex>
#include <vector>
#include <numeric>
#include <cstdlib>
#include <cstdio>

#define NTRIAL 11
double test(int N) {
	double err, t0 = 0, t1 = 0;
	for (int n = 0; n < NTRIAL + 1; n++) {
		std::vector<std::complex<double>> Z1(N), Z2(N);
		for (int n = 0; n < N; n++) {
			Z1[n] = Z2[n] = std::complex<double>(rand1(), rand1());
		}
		Z1[3] = Z2[3] = 1.0;
		t0 += fftw(Z2.data(), N);
		t1 += fft_1d(Z1.data(), N);
		if (n == 0) {
			t1 = t0 = 0.0;
		}
		err = 0.0;
		for (int n = 0; n < N; n++) {
			err += std::abs(Z1[n] - Z2[n]);
			//	printf("%3i | %15e %15e  | %15e %15e  | %15e %15e |\n", n, Z1[n].real(), Z1[n].imag(), Z2[n].real(),
			//			Z2[n].imag(), Z2[n].real() - Z1[n].real(), Z2[n].imag() - Z1[n].imag());
		}
		err /= N;
		//	abort();
	}
	printf("%20i %20e %20e %20e %20f\n", N, err, t0, t1, t0 / (1.0e-10 + t1));
	return t1;
}

void * operator new(std::size_t n) {
	void* memptr;
	posix_memalign(&memptr, 32, round_up(n, 32));
	return memptr;
}

void operator delete(void * p) {
	free(p);
}

void *operator new[](std::size_t n) {
	void* memptr;
	posix_memalign(&memptr, 32, round_up(n, 32));
	return memptr;
}

void operator delete[](void *p) {
	free(p);
}

bool usenum(int N) {
	int M = N;
	int cnt = 0;
	while (M > 1) {
		if (M % 5 == 0) {
			M /= 5;
		} else if (M % 3 == 0) {
			M /= 3;
		} else if (M % 2 == 0) {
			M /= 2;
		} else {
			return false;
		}
		cnt++;
	}
	return cnt > 1;
}

#include <fenv.h>

int main() {
	//feenableexcept(FE_DIVBYZERO);
	//feenableexcept(FE_INVALID);
//	feenableexcept(FE_OVERFLOW);
	int N = 81;
	int j = 0;
	for (int i = 0; i < N; i++) {
	if (i != N - 1) {
			int k = N / 3;
			j++;
			while (2 * k < j) {
				j -= 2 * k;
				k /= 3;
			}
			j += k - 1;
		}
	}
//	return 0;

	std::vector<int> nums;
	double tm = 0.0;
	int R = 13;
	for (int N = R * R; N <= 100000000; N *= R) {
		tm += test(N);

	}
	return 0;
}
