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
#include <fft86/scramble.hpp>
#include <fft86/twiddles.hpp>

#include <algorithm>
#include <complex>
#include <vector>
#include <numeric>
#include <cstdlib>
#include <cstdio>

void fft_pfa(std::complex<double>* X, int N1, int N2);

void fft_complex(double* X, double* Y, int N);

#define NTRIAL 101
void pfa(std::complex<double>* X, int N);

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
			//		printf("%3i | %15e %15e  | %15e %15e  | %15e %15e |\n", n, Z1[n].real(), Z1[n].imag(), Z2[n].real(),
			//				Z2[n].imag(), Z2[n].real() - Z1[n].real(), Z2[n].imag() - Z1[n].imag());
		}
		err /= N;
		//	abort();
	}
	printf("%20i %20e %20e %20e %20f - ", N, err, t0, t1, t0 / (1.0e-10 + t1));
	int M = N;
	int n = 2;
	while (M > 1) {
		for (; n <= M; n++) {
			if (M % n == 0) {
				printf("%i", n);
				int cnt = 0;
				while (M % n == 0) {
					M /= n;
					cnt++;
				}
				if (cnt > 1) {
					printf("^%i", cnt);
				}
				if (M != 1) {
					printf("*");
				}
			}
		}
	}
	printf("\n");
	return t1;
}

void test_pfa(int N) {
	double err, t0 = 0, t1 = 0;
	for (int n = 0; n < NTRIAL + 1; n++) {
		std::vector<std::complex<double>> Z1(N), Z2(N);
		for (int n = 0; n < N; n++) {
			Z1[n] = Z2[n] = std::complex<double>(0,0);
		}
		Z1[7] = Z2[7] = 1.0;
		t0 += fftw(Z2.data(), N);
		timer tm;
		tm.start();
		pfa(Z1.data(), N);
		tm.stop();
		t1 += tm.read();
		if (n == 0) {
			t1 = t0 = 0.0;
		}
		err = 0.0;
		for (int n = 0; n < N; n++) {
			err += std::abs(Z1[n] - Z2[n]);
			printf("%3i | %15e %15e  | %15e %15e  | %15e %15e |\n", n, Z1[n].real(), Z1[n].imag(), Z2[n].real(),
					Z2[n].imag(), Z2[n].real() - Z1[n].real(), Z2[n].imag() - Z1[n].imag());
		}
		err /= N;
	//	abort();
	}

	printf("%20i %20e %20e %20e %20f - ", N, err, t0, t1, t0 / (1.0e-10 + t1));
	printf("\n");
}

#include <sfft.hpp>

int main() {
	test_pfa(2*2*3);
//feenableexcept(FE_DIVBYZERO);
//feenableexcept(FE_INVALID);
//	feenableexcept(FE_OVERFLOW);
	return 0;
}
