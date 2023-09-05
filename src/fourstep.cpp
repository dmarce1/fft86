#include <fft86/scramble.hpp>
#include <fft86/timer.hpp>
#include <fft86/twiddles.hpp>
#include <fft86/util.hpp>
#include <fft86/vec.hpp>
#include <fft86/defs.hpp>

#include <array>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <complex>
#include <cstring>
#include <vector>
#include <sfft.hpp>

#define L1SIZE (32 * 1024)
#define L2SIZE (512 * 1024)
#define L3SIZE (16384 * 1024)

struct params_t {
	int N;
	const std::vector<double>* WR;
	const std::vector<double>* WI;
	const double* wr;
	const double* wi;
	const int* P;
};

template<int N1>
void fft_1d_batch_dit(const params_t* params, double* Z, int N, int NLO, int leve);

void perf_shuf(double* Z, int N2, int NLO) {
	const int NLO4 = round_down(NLO, v4df::size());
	const int NLO2 = round_down(NLO, v2df::size());
	const int& NLO1 = NLO;
	for (int k2 = 0; k2 < N2; k2++) {
		for (int ilo = 0; ilo < NLO4; ilo += v4df::size()) {
			auto* r = Z + 2 * (ilo + NLO * k2);
			auto* i = r + 4;
			sfft_perf_shuf_w4(r, i);
		}
		if (NLO2 > NLO4) {
			auto* r = Z + 2 * (NLO4 + NLO * k2);
			auto* i = r + 2;
			sfft_perf_shuf_w2(r, i);
		}
	}
}

void inv_perf_shuf(double* Z, int N2, int NLO) {
	const int NLO4 = round_down(NLO, v4df::size());
	const int NLO2 = round_down(NLO, v2df::size());
	const int& NLO1 = NLO;
	for (int k2 = 0; k2 < N2; k2++) {
		for (int ilo = 0; ilo < NLO4; ilo += v4df::size()) {
			auto* r = Z + 2 * (ilo + NLO * k2);
			auto* i = r + 4;
			sfft_inv_perf_shuf_w4(r, i);
		}
		if (NLO2 > NLO4) {
			auto* r = Z + 2 * (NLO4 + NLO * k2);
			auto* i = r + 2;
			sfft_inv_perf_shuf_w2(r, i);
		}
	}
}
/*
 template<int N1>
 void fft_1d_batch_dif(double* Z, int N, int NLO, int level) {
 const auto* W = get_twiddles(N);
 const int NLO4 = round_down(NLO, v4df::size());
 const int NLO2 = round_down(NLO, v2df::size());
 const int& NLO1 = NLO;
 const int N2 = N / N1;
 const int s = 2 * NLO * N2;
 for (int k2 = 0; k2 < N2; k2++) {
 if (level == 0) {
 for (int n1 = 0; n1 < N1; n1++) {
 inv_perf_shuf(Z + 2 * NLO * k2 + s * n1, 1, NLO);
 }
 }
 std::array<double, 8 * N1> w;
 for (int n = 0; n < N1; n++) {
 for (int l = 0; l < 4; l++) {
 w[8 * n + l] = W[n * k2].real();
 w[8 * n + l + 4] = W[n * k2].imag();
 }
 }
 for (int ilo = 0; ilo < NLO4; ilo += v4df::size()) {
 double* x = Z + 2 * (ilo + NLO * k2);
 double* y = x + 4;
 sfft_complex_dif_w4(x, y, s, N1, w.data());
 }
 if (NLO2 > NLO4) {
 double* x = Z + 2 * (NLO4 + NLO * k2);
 double* y = x + 2;
 for (int n = 0; n < N1; n++) {
 for (int l = 0; l < 2; l++) {
 w[4 * n + l] = W[n * k2].real();
 w[4 * n + l + 2] = W[n * k2].imag();
 }
 }
 sfft_complex_dif_w2(x, y, s, N1, w.data());
 }
 if (NLO1 > NLO4 && NLO1 != NLO2) {
 double* x = Z + 2 * (NLO2 + NLO * k2);
 double* y = x + 1;
 for (int n = 0; n < N1; n++) {
 w[2 * n] = W[n * k2].real();
 w[2 * n + 1] = W[n * k2].imag();
 }
 sfft_complex_dif_w1(x, y, s, N1, w.data());
 }
 }
 if (N2 > 1) {
 for (int n1 = 0; n1 < N1; n1++) {
 fft_1d_batch_dif<N1>(Z + 2 * n1 * NLO * N2, N2, NLO, level + 1);
 }
 }
 }*/

template<int N1>
void fft_1d_batch_dit(double* Z, int N, int NLO, int level) {
	const auto* W = get_twiddles(N);
	const int NLO4 = round_down(NLO, v4df::size());
	const int NLO2 = round_down(NLO, v2df::size());
	const int& NLO1 = NLO;
	const int N2 = N / N1;
	if (N2 > 1) {
		for (int n1 = 0; n1 < N1; n1++) {
			fft_1d_batch_dit<N1>(Z + 2 * n1 * NLO * N2, N2, NLO, level + 1);
		}
	}
	const int s = 2 * NLO * N2;
	for (int k2 = 0; k2 < N2; k2++) {
		if (N2 == 1) {
			for (int n1 = 0; n1 < N1; n1++) {
				inv_perf_shuf(Z + 2 * NLO * k2 + s * n1, 1, NLO);
			}
		}
		std::array<double, 8 * N1> w;
		for (int n = 0; n < N1; n++) {
			for (int l = 0; l < 4; l++) {
				w[8 * n + l] = W[n * k2].real();
				w[8 * n + l + 4] = W[n * k2].imag();
			}
		}
		for (int ilo = 0; ilo < NLO4; ilo += v4df::size()) {
			double* x = Z + 2 * (ilo + NLO * k2);
			double* y = x + 4;
			sfft_complex_dit_w4(x, y, s, N1, w.data());
		}
		if (NLO2 > NLO4) {
			double* x = Z + 2 * (NLO4 + NLO * k2);
			double* y = x + 2;
			for (int n = 0; n < N1; n++) {
				for (int l = 0; l < 2; l++) {
					w[4 * n + l] = W[n * k2].real();
					w[4 * n + l + 2] = W[n * k2].imag();
				}
			}
			sfft_complex_dit_w2(x, y, s, N1, w.data());
		}
		if (NLO1 > NLO4 && NLO1 != NLO2) {
			double* x = Z + 2 * (NLO2 + NLO * k2);
			double* y = x + 1;
			for (int n = 0; n < N1; n++) {
				w[2 * n] = W[n * k2].real();
				w[2 * n + 1] = W[n * k2].imag();
			}
			sfft_complex_dit_w1(x, y, s, N1, w.data());
		}
		if (level == 0) {
			for (int n1 = 0; n1 < N1; n1++) {
				perf_shuf(Z + 2 * NLO * k2 + s * n1, 1, NLO);
			}
		}
	}
}
/*
template<int R>
void apply_twiddles(std::complex<double>* Z, int N1, int N2) {
	const auto& P = get_digit_reversal(R, N2);
	const auto& Wr = *get_6step_cos_twiddles(N2, N1);
	const auto& Wi = *get_6step_sin_twiddles(N2, N1);
	const int N = N1 * N2;
	const int end4 = round_down(N1, v4df::size());
	const int end2 = round_down(N1, v2df::size());
	for (int k2 = 0; k2 < N2; k2++) {
		const int k2r = P[k2];
		const int k20 = k2;
		for (int n1 = 0; n1 < end4; n1 += v4df::size()) {
			double* X = (double*) Z + 2 * (n1 + N1 * k2r);
			double* Y = (double*) X + v4df::size();
			v4df x, y, tmp, wr, wi;
			wr.load(Wr[k20].data() + n1);
			wi.load(Wi[k20].data() + n1);
			x.load(X);
			y.load(Y);
			tmp = x;
			x = tmp * wr - y * wi;
			y = tmp * wi + y * wr;
			x.store(X);
			y.store(Y);
			sfft_perf_shuf_w4(X, Y);
		}
		if (end2 > end4) {
			const int& n1 = end4;
			double* X = (double*) Z + 2 * (n1 + N1 * k2r);
			double* Y = (double*) X + v2df::size();
			v2df x, y, tmp, wr, wi;
			wr.load(Wr[k20].data() + n1);
			wi.load(Wi[k20].data() + n1);
			x.load(X);
			y.load(Y);
			tmp = x;
			x = tmp * wr - y * wi;
			y = tmp * wi + y * wr;
			x.store(X);
			y.store(Y);
			sfft_perf_shuf_w2(X, Y);
		}
		if (N1 > end4 && N1 != end2) {
			const int& n1 = end2;
			double* X = (double*) Z + 2 * (n1 + N1 * k2r);
			double* Y = (double*) X + v1df::size();
			double x, y, tmp, wr, wi;
			wr = *(Wr[k20].data() + n1);
			wi = *(Wi[k20].data() + n1);
			x = *X;
			y = *Y;
			tmp = x;
			x = tmp * wr - y * wi;
			y = tmp * wi + y * wr;
			*X = x;
			*Y = y;
		}
	}
}*/

extern "C" void fft_1d_batch_dif(double* Z, const double* W, int N1, int N, int NLO);
extern "C" void fft_1d_batch_dit(double* Z, const double* W, int N1, int N, int NLO);
extern "C" void apply_twiddles(double* Z, const double* W, const double*, int R, int N1, int N2);
timer tm1, tm2, tm3, tm4;

template<int R>
void fft_1d(std::complex<double>* Z, int N) {
	int ilogb = std::lround(std::log(N) / std::log(R)) >> 1;
	int N1 = 1;
	for (int i = 0; i < ilogb; i++) {
		N1 *= R;
	}
	int N2 = N / N1;

	tm2.start();
	const auto* W1 = get_twiddles(N1);
	fft_1d_batch_dif((double*) Z, (const double*) W1, R, N1, N2);
	tm2.stop();

	const auto& Wr = *get_6step_cos_twiddles(N2, N1);
	const auto& Wi = *get_6step_sin_twiddles(N2, N1);
	tm3.start();
	apply_twiddles((double*) Z, Wr.data(), Wi.data(), R, N2, N1);
	tm3.stop();

	tm4.start();
	scramble(Z, R, N);
	tm4.stop();

	tm2.start();
	const auto* W2 = get_twiddles(N2);
	fft_1d_batch_dit((double*) Z, (const double*) W2, R, N2, N1);
	tm2.stop();

	double tinv = 1.0 / (tm1.read() + tm2.read() + tm3.read() + tm4.read());
//	printf( "%e %e %e %e\n", tm1.read() * tinv, tm2.read() * tinv, tm3.read() * tinv, tm4.read() * tinv);

}

double fft_1d(std::complex<double>* Z, int N) {
	timer tm;
	tm.start();
	int R = 4;
	while (N % R != 0) {
		R++;
	}
	switch (R) {
	case 2:
		fft_1d<2>(Z, N);
		break;
	case 3:
		fft_1d<3>(Z, N);
		break;
	case 4:
		fft_1d<4>(Z, N);
		break;
	case 5:
		fft_1d<5>(Z, N);
		break;
	case 7:
		fft_1d<7>(Z, N);
		break;
	case 8:
		fft_1d<8>(Z, N);
		break;
	case 9:
		fft_1d<9>(Z, N);
		break;
	case 11:
		fft_1d<11>(Z, N);
		break;
	case 13:
		fft_1d<13>(Z, N);
		break;
	case 25:
		fft_1d<25>(Z, N);
		break;
	case 32:
		fft_1d<32>(Z, N);
		break;
	default:
		assert(false);
		abort();
	}
	tm.stop();
	return tm.read();
}
