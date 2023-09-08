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
#include <fft86/timer.hpp>

#include <algorithm>
#include <complex>
#include <vector>
#include <numeric>
#include <cstdlib>
#include <cstdio>

void fft_pfa(std::complex<double>* X, int N1, int N2);

void fft_complex(double* X, double* Y, int N);

#define NTRIAL 101
double test(int N) {
	double err, t0 = 0, t1 = 0;
	for (int n = 0; n < NTRIAL + 1; n++) {
		std::vector<std::complex<double>> Z1(N), Z2(N);
		std::vector<double> X(N), Y(N);
		for (int n = 0; n < N; n++) {
			Z1[n] = Z2[n] = std::complex<double>(rand1(), rand1());
		}
		Z1[1] = Z2[1] = 1.0;
		for (int n = 0; n < N; n++) {
			X[n] = Z1[n].real();
			Y[n] = Z1[n].imag();
		}
		t0 += fftw(Z2.data(), N);
		timer tm;
		tm.start();
		fft_complex(X.data(), Y.data(), N);
		tm.stop();
		t1 += tm.read();
		for (int n = 0; n < N; n++) {
			Z1[n].real(X[n]);
			Z1[n].imag(Y[n]);
		}
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

void fft_prime_factor_algorithm(double* X, const std::vector<int>& N);

double test_pfa(std::vector<int> Ns) {
	int N = std::reduce(Ns.begin(), Ns.end(), 1, std::multiplies<int>());
	double err, t0 = 0, t1 = 0;
	for (int n = 0; n < NTRIAL + 1; n++) {
		std::vector<std::complex<double>> Z1(N), Z2(N);
		for (int n = 0; n < N; n++) {
			Z1[n] = Z2[n] = std::complex<double>(rand1(), rand1());
		}
		Z1[0] = Z2[0] = 1.0;
		t0 += fftw(Z2.data(), N);
		timer tm;
		tm.start();
		fft_prime_factor_algorithm((double*) Z1.data(), Ns);
		tm.stop();
		t1 += tm.read();
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
		//abort();
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

void test2(int N1, int N2) {
	int n1 = 1;
	std::vector<bool> visited(N1, false);
	for (int n1 = 0; n1 < N1; n1++) {
		if (!visited[n1]) {
			printf("\n");
			int n = n1;
			do {
				int next = (n * N2) % N1;
				printf("%i %i\n", n, next);
				n = next;
				visited[n] = true;
			} while (n != n1);
		}
	}
}

std::vector<std::pair<int, int>> compute_factors(int N);

void sort_nonsquare(double* X, int* N, int nN);
void sort(double* X, double* Y, int N, int NLO = 1);

int main() {
	test(2*2*14*13*17*19);
	//feenableexcept(FE_DIVBYZERO);
	//feenableexcept(FE_INVALID);
//	feenableexcept(FE_OVERFLOW);
	return 0;
}
