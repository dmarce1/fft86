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

void dct3(double* Y, int N) {
	std::vector<std::complex<double>> Z(N / 2 + 1);
	std::vector<double> X(N);
	Z[0] = Y[0];
	Z[N / 2] = sqrt(0.5)*Y[N / 2];
	for (int n = 1; n < (N+1) / 2; n++) {
		Z[n] = 0.5 * (Y[n] - std::complex<double>(0, 1) * Y[N - n]) * std::polar(1.0, M_PI * 0.5 * n / N);
	}
	fftw_inv_real(Z, X);
	for (int n = 0; n <= N / 2; n++) {
		Y[2 * n] = X[n];
		Y[2 * n + 1] = X[N - n - 1];
	}
}

void dst3(double* X, int N) {
	std::vector<double> Y(N);
	for (int k = 1; k < N - k; k++) {
		std::swap(X[k], X[N - k]);
	}
	dct3(X, N);
	for (int k = 1; k < N; k += 2) {
		X[k] = -X[k];
	}
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

void fft_skew2(double* X, int N) {
	if (N % 2 == 0) {
		std::vector<double> xp(N);
		std::vector<double> xm(N);
		for (int n = 0; n < N / 2; n++) {
			xp[n] = n == 0 ? X[0] : (X[n] - X[(N - n) % N]);
			xm[n] = n == 0 ? 0.0 : -(X[n] + X[(N - n) % N]);
		}
		dct3(xp.data(), N / 2);
		dst3(xm.data(), N / 2);
		const double xny = X[N / 2];
		for (int n = 0; n < N / 2; n++) {
			X[n] = xp[n];
			X[N - n - 1] = xm[n] - std::pow(-1, n) * xny;
		}
	} else {
		std::vector<double> y(N);
		const int dn = N / 2;
		for (int n = 0; n < N; n++) {
			y[n] = X[n];
			y[n] *= std::pow(-1, n);
		}
		std::vector<std::complex<double>> z((N + 1) / 2);
		fftw_real(y, z);
		std::complex<double> z0(0.0, 0.0);
		for (int n = 0; n < N; n++) {
			z0 += X[n] * std::polar(1.0, -M_PI * n / N);
		}
		for (int n = 0; n < N / 2; n++) {
			X[N / 2 - n] = z[(n + dn) % (N / 2)].real();
			X[N / 2 + n] = -z[(n + dn) % (N / 2)].imag();
		}
		X[0] = z0.real();
		X[N - 1] = z0.imag();
		X[N / 2] = z[0].real();
	}
}

int main() {
	int N = 16;
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
