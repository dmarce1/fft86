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

void fft_prime_factor_algorithm(double* X, const std::vector<int>& N);

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

#include <fftw3.h>
#include <unordered_map>

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

double fftw_dctIII(double* x, int N) {
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, double*> in;
	static std::unordered_map<int, double*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = out[N] = (double*) fftw_malloc(sizeof(double) * N);
		plans[N] = fftw_plan_r2r_1d(N, in[N], out[N], FFTW_REDFT01, FFTW_MEASURE);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N; n++) {
		i[n] = x[n];
	}
	timer tm;
	tm.start();
	fftw_execute(plans[N]);
	tm.stop();
	for (int n = 0; n < N; n++) {
		x[n] = o[n];
	}
	return tm.read();
}

double fftw_dctII(double* x, int N) {
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, double*> in;
	static std::unordered_map<int, double*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = out[N] = (double*) fftw_malloc(sizeof(double) * N);
		plans[N] = fftw_plan_r2r_1d(N, in[N], out[N], FFTW_REDFT10, FFTW_MEASURE);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N; n++) {
		i[n] = x[n];
	}
	timer tm;
	tm.start();
	fftw_execute(plans[N]);
	tm.stop();
	for (int n = 0; n < N; n++) {
		x[n] = o[n];
	}
	return tm.read();
}

double fftw_dstIII(double* x, int N) {
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, double*> in;
	static std::unordered_map<int, double*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = out[N] = (double*) fftw_malloc(sizeof(double) * N);
		plans[N] = fftw_plan_r2r_1d(N, in[N], out[N], FFTW_RODFT01, FFTW_MEASURE);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N; n++) {
		i[n] = x[n];
	}
	timer tm;
	tm.start();
	fftw_execute(plans[N]);
	tm.stop();
	for (int n = 0; n < N; n++) {
		x[n] = o[n];
	}
	return tm.read();
}

double fftw_real(const std::vector<double>& xin, std::vector<std::complex<double>>& xout) {
	const int N = xin.size();
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, double*> in;
	static std::unordered_map<int, fftw_complex*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = (double*) malloc(sizeof(double) * N);
		out[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * (N / 2 + 1));
		plans[N] = fftw_plan_dft_r2c_1d(N, in[N], out[N], FFTW_MEASURE);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N; n++) {
		i[n] = xin[n];
	}
	timer tm;
	tm.start();
	fftw_execute(plans[N]);
	tm.stop();
	for (int n = 0; n < N / 2 + 1; n++) {
		xout[n].real((o[n][0]));
		xout[n].imag((o[n][1]));
	}
	return tm.read();
}

double fftw_inv_real(const std::vector<std::complex<double>>& xin, std::vector<double>& xout) {
	const int N = xout.size();
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, double*> out;
	static std::unordered_map<int, fftw_complex*> in;
	if (plans.find(N) == plans.end()) {
		out[N] = (double*) malloc(sizeof(double) * N);
		in[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * (N / 2 + 1));
		plans[N] = fftw_plan_dft_c2r_1d(N, in[N], out[N], FFTW_MEASURE);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N / 2 + 1; n++) {
		i[n][0] = ((xin[n].real()));
		i[n][1] = ((xin[n].imag()));
	}
	timer tm;
	tm.start();
	fftw_execute(plans[N]);
	tm.stop();
	for (int n = 0; n < N; n++) {
		xout[n] = o[n];
	}
	return tm.read();
}

void dst3(double* Y, int N) {

	std::vector<std::complex<double>> Z(N / 2 + 1);
	std::vector<double> X(N);
	X[0] = Y[0];
	X[N / 2] = sqrt(0.5) * Y[N / 2];
	for (int n = 1; n < N - n; n++) {
		const double c = cos(0.5 * M_PI * n / N);
		const double s = sin(0.5 * M_PI * n / N);
		X[N - n] = 0.5 * (Y[N - n] * (c - s) + Y[n] * (s + c));
		X[n] = 0.5 * (Y[N - n] * (c + s) + Y[n] * (s - c));
	}
	fftw_real(X, Z);
	Y[0] = Z[0].real();
	Y[N / 2] = Z[N / 2].real();
	for (int n = 1; n <= N / 2; n++) {
		Y[2 * n] = Z[n].real() + Z[n].imag();
		Y[2 * n - 1] = Z[n].imag() - Z[n].real();
	}
}

void dct3(double* Y, int N) {
	std::vector<std::complex<double>> Z(N / 2 + 1);
	std::vector<double> X(N);
	X[0] = Y[0];
	X[N / 2] = sqrt(0.5) * Y[N / 2];
	for (int n = 1; n < N - n; n++) {
		const double c = cos(0.5 * M_PI * n / N);
		const double s = sin(0.5 * M_PI * n / N);
		X[n] = 0.5 * (Y[n] * (c + s) + Y[N - n] * (s - c));
		X[N - n] = 0.5 * (Y[n] * (c - s) + Y[N - n] * (s + c));
	}
	fftw_real(X, Z);
	Y[0] = Z[0].real();
	Y[N / 2] = Z[N / 2].real();
	for (int n = 1; n <= N / 2; n++) {
		Y[2 * n] = Z[n].real() + Z[n].imag();
		Y[2 * n - 1] = Z[n].real() - Z[n].imag();
	}
}

#include <cstring>

void fft_skew2(double* x, int N) {
	std::vector<double> X(x, x + N);
	const auto* W2 = get_twiddles(2 * N);
	if (N % 2 == 0) {
		auto P = get_digit_reversal(2, N);
		for (int n = 0; n < N; n++) {
			if (P[n] > n) {
				std::swap(X[P[n]], X[n]);
			}
		}
		std::vector<std::complex<double>> Z(N / 2 + 1);
		double T1 = 0.0, T2 = 0.0;
		for (int n = 0; n < N / 2; n++) {
			T1 += X[n] * W2[P[n]].real();
			T2 += X[n + N / 2] * W2[P[n + N / 2]].real();
		}
		for (int n = 0; n < N; n++) {
			X[n] *= -2.0 * W2[P[n]].imag();
		}
		for (int n = 0; n < N; n++) {
			if (P[n] > n) {
				std::swap(X[P[n]], X[n]);
			}
		}
		fftw_real(X, Z);
		for (int n = 1; n < N / 2; n++) {
			X[n] = Z[n].real();
			X[N - n] = Z[n].imag();
		}
		X[0] = Z[0].real();
		X[N / 2] = Z[N / 2].real();
		double T3 = X[N - 1];
		X[N - 1] = -0.5 * X[0];
		X[0] = T1 + T2;
		double T4;
		for (int k1 = 1; k1 < N / 2 - 1; k1++) {
			T4 = X[N - k1 - 1];
			X[N - k1 - 1] = -X[k1] + X[N - k1];
			X[k1] = T3 + X[k1 - 1];
			T3 = T4;
		}
		X[N / 2] = 0.5 * X[N / 2];
		X[N / 2 - 1] = T1 - T2;
	} else {
		auto P = get_digit_reversal(5, N);
		for (int n = 0; n < N; n++) {
			if (P[n] > n) {
				std::swap(X[P[n]], X[n]);
			}
		}
		std::vector<double> y(N);
		const int dn = N / 2;
		for (int n = 0; n < N; n++) {
			y[n] = X[n];
			y[n] *= std::pow(-1, P[n]);
		}
		std::vector<std::complex<double>> z((N + 1) / 2);
		for (int n = 0; n < N; n++) {
			if (P[n] > n) {
				std::swap(y[P[n]], y[n]);
			}
		}
		fftw_real(y, z);
		std::complex<double> z0(0.0, 0.0);
		for (int n = 0; n < N; n++) {
			z0 += X[n] * W2[P[n]];
		}
		X[0] = z[0].real();
		for (int n = 1; n < N - n; n++) {
			X[n] = z[n].real();
			X[N - n] = z[n].imag();
		}
		for (int n = 1; n <= N / 4; n++) {
			std::swap(X[N - n], X[N / 2 + n]);
			X[N - n] *= -1;
			X[N / 2 + n] *= -1;
			std::swap(X[n], X[N / 2 - n]);
		}
		X[0] = z0.real();
		X[N - 1] = z0.imag();
		X[N / 2] = z[0].real();
	}
	std::memcpy(x, X.data(), sizeof(double) * N);
}

void fft_skew1(double* X, int N) {
	std::vector<std::complex<double>> Y(N);
	for (int k = 0; k < N; k++) {
		Y[k] = 0.0;
		for (int n = 0; n < N; n++) {
			const auto w = std::polar(1.0, -2.0 * M_PI * n * (k + 0.5) / N);
			Y[k] += X[n] * w;
		}
	}
	for (int k = 0; k < N - k; k++) {
		X[k] = Y[k].real();
		if (N - k - 1 != k) {
			X[N - k - 1] = Y[k].imag();
		}
	}
}

void dft_skew(double* X, int N) {
	std::vector<std::complex<double>> Y(N);
	for (int k = 0; k < N; k++) {
		Y[k] = 0.0;
		for (int n = 0; n < N; n++) {
			const auto w = std::polar(1.0, -2.0 * M_PI * n * (k + 0.5) / N);
			Y[k] += X[n] * w;
		}
	}
	for (int k = 0; k < N - k; k++) {
		X[k] = Y[k].real();
		if (N - k - 1 != k) {
			X[N - k - 1] = Y[k].imag();
		}
	}
}

void dft(double* X, int N) {
	std::vector<std::complex<double>> Y(N);
	for (int k = 0; k < N; k++) {
		Y[k] = 0.0;
		for (int n = 0; n < N; n++) {
			const auto w = std::polar(1.0, -2.0 * M_PI * n * k / N);
			Y[k] += X[n] * w;
		}
	}
	for (int k = 0; k < N - k; k++) {
		X[k] = Y[k].real();
		if (N - k - 1 != k) {
			X[N - k] = Y[k].imag();
		}
	}
}

#include <sfft.hpp>

int main() {
	int N = 25;
	std::vector<double> X(N);
	std::vector<double> Y(N);
	for (int n = 0; n < N; n++) {
		X[n] = rand1();
	}
	X[2] = 1.0;
	for (int n = 0; n < N; n++) {
		Y[n] = X[n];
	}
	auto P = get_digit_reversal(3, N);
	for (int n = 0; n < N; n++) {
		if (P[n] > n) {
			//	std::swap(X[P[n]], X[n]);
		}
	}

	dft_skew(Y.data(), N);
	fft_skew2(X.data(), N);
	for (int n = 0; n < N; n++) {
		printf("%i %e %e %e\n", n, X[n], Y[n], X[n] - Y[n]);
	}
	/*
	 for (int N = 65536; N < 1024 * 1024 * 1024; N++) {
	 int M = N;
	 int nfac = 0;
	 int facmax;
	 for (int n = 4; n <= N; n++) {
	 while (M % n == 0) {
	 facmax = n;
	 M /= n;
	 nfac++;
	 }
	 }
	 if (facmax <= 32 && nfac >= 2) {
	 test(N);
	 N *= 11;
	 N /= 10;
	 N--;
	 }
	 }*/
//feenableexcept(FE_DIVBYZERO);
//feenableexcept(FE_INVALID);
//	feenableexcept(FE_OVERFLOW);
	return 0;
}
