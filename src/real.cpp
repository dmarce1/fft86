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

enum sort_type {
	TOP, MIDDLE
};

int reverse_digits(int* digs, int ndigs, int index) {
	int k = index;
	int j = 0;
	for (int l = 0; l < ndigs; l++) {
		std::div_t q = std::div(k, digs[l]);
		j *= digs[l];
		j += q.rem;
		k = q.quot;
	}
	return j;
}

void reverse_digits(double* X, int* digs, int ndigs, int M) {
	const int N = std::reduce(digs, digs + ndigs, 1, std::multiplies<int>());
	static thread_local std::vector<bool> visited;
	visited.resize(N);
	std::fill(visited.begin(), visited.end(), false);
	for (int i = 1; i < N - 1; i++) {
		if (!visited[i]) {
			int j = i;
			do {
				j = reverse_digits(digs, ndigs, j);
				for (int m = 0; m < M; m++) {
					std::swap(X[i * M + m], X[j * M + m]);
				}
				visited[j] = true;
			} while (j != i);
		}
	}
}

void reverse_digits(double* X, int* digs, int ndigs, int beg, int end) {
	const int NHI = std::reduce(digs, digs + beg, 1, std::multiplies<int>());
	const int NMID = std::reduce(digs + beg, digs + end, 1, std::multiplies<int>());
	const int NLO = std::reduce(digs + end, digs + ndigs, 1, std::multiplies<int>());
	for (int i = 0; i < NHI; i++) {
		reverse_digits(X + i * NLO * NMID, digs + beg, end - beg, NLO);
	}
}

void sfft_real(double* X, int s, int N1, int width) {
	switch (width) {
	case 4:
		return sfft_real_w4(X, s, N1);
	case 2:
		return sfft_real_w2(X, s, N1);
	case 1:
		return sfft_real_w1(X, s, N1);
	}
}

void sfft_skew(double* X, int s, int N1, int width) {
	switch (width) {
	case 4:
		return sfft_skew_w4(X, s, N1);
	case 2:
		return sfft_skew_w2(X, s, N1);
	case 1:
		return sfft_skew_w1(X, s, N1);
	}
}

void sfft_real_complex(double* X, double* Y, int s, int N1, double* w, int width) {
	switch (width) {
	case 4:
		return sfft_real_complex_w4(X, Y, s, N1, w);
	case 2:
		return sfft_real_complex_w2(X, Y, s, N1, w);
	case 1:
		return sfft_real_complex_w1(X, Y, s, N1, w);
	}
}

void sfft_complex(double* X, double* Y, int s, int N1, double* w, int width) {
	switch (width) {
	case 4:
		return sfft_real_complex_w4(X, Y, s, N1, w);
	case 2:
		return sfft_real_complex_w2(X, Y, s, N1, w);
	case 1:
		return sfft_real_complex_w1(X, Y, s, N1, w);
	}
}

void fft_real(double* X, int N, int* digb, int* dige, int M) {
	const int N1 = *digb;
	const int N2 = N / N1;
	if (N2 > 1) {
		for (int n1 = 0; n1 < N1; n1++) {
			fft_real(X + n1 * N2, N2, digb + 1, dige, M);
		}
	}
	for (int w = 4, m = 0; w >= 1; w >>= 1) {
		for (; m <= M - w; m += w) {
			sfft_real(X + m, M * N2, N1, w);
		}
	}
	if (N2 / 2) {
		for (int w = 4, m = 0; w >= 1; w >>= 1) {
			for (; m <= M - w; m += w) {
				sfft_skew(X + M * N2 / 2 + m, M * N2, N1, w);
			}
		}
	}
	if (N2 / 2 > 1) {
		static thread_local std::vector<double> wn;
		wn.resize(2 * N1);
		const auto* WN = get_twiddles(N);
		for (int k2 = 1; k2 < N2 / 2; k2++) {
			for (int w = 4, m = 0; w >= 1; w >>= 1) {
				for (int n1 = 0; n1 < N1; n1++) {
					for (int l = 0; l < w; l++) {
						wn[2 * w * n1 + l] = WN[n1 * k2].real();
						wn[2 * w * n1 + l + w] = WN[n1 * k2].imag();
					}
				}
				for (; m <= M - w; m += w) {
					sfft_real_complex(X + M * k2 + m, X + M * (N2 - k2) + m, M * N2, N1, wn.data(), w);
				}
			}
		}
	}
}

void fft_post_skew(double* X, int N, int M, std::pair<double, double> in) {
	const auto* W2 = get_twiddles(2 * N);
	if (N % 2 == 0) {
		double T1 = in.first;
		double T2 = in.second;
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
		std::complex<double> z0(in.first, in.second);
		for (int n = 1; n <= N / 4; n++) {
			std::swap(X[N - n], X[N / 2 + n]);
			X[N - n] = -X[N - n];
			X[N / 2 + n] = -X[N / 2 + n];
			std::swap(X[n], X[N / 2 - n]);
		}
		X[N / 2] = X[0];
		X[0] = z0.real();
		X[N - 1] = z0.imag();
	}
}

std::pair<double, double> fft_pre_skew(double* X, int N, int* digb, int* dige, int M) {
	const auto* W2 = get_twiddles(2 * N);
	if (N % 2 == 0) {
		std::vector<std::complex<double>> Z(N / 2 + 1);
		double T1 = 0.0, T2 = 0.0;
		for (int n = 0; n < N / 2; n++) {
			const auto pn = reverse_digits(digb, dige - digb, n);
			T1 += X[n] * W2[pn].real();
			X[n] *= -2.0 * W2[pn].imag();
		}
		for (int n = N / 2; n < N; n++) {
			const auto pn = reverse_digits(digb, dige - digb, n);
			T2 += X[pn] * W2[pn].real();
			X[n] *= -2.0 * W2[pn].imag();
		}
		return std::make_pair(T1, T2);
	} else {
		std::complex<double> z0(0.0, 0.0);
		for (int n = 0; n < N; n++) {
			const auto pn = reverse_digits(digb, dige - digb, n);
			z0 += X[n] * W2[pn];
			X[n] *= std::pow(-1, pn);
		}
		return std::make_pair(z0.real(), z0.imag());
	}
}

void form_complex(double* X, int N1, int N2, int M) {
	for (int k1 = 0; k1 < N2; k1++) {
		for (int k2 = 1; k2 < N - k2; k2++) {
			const auto lri = X + N2 * M * k1 + k2 * M;
			const auto rri = X + N2 * M * k1 + (N2 - k2) * M;
			const auto lii = X + N2 * M * (N1 - k1) + k2 * M;
			const auto rii = X + N2 * M * (N1 - k1) + (N2 - k2) * M;
			for (int m = 0; m < M; m++) {
				const auto x1 = X[lri + m] - X[rii + m];
				const auto x2 = X[lri + m] + X[rii + m];
				const auto y1 = X[lii + m] + X[rri + m];
				const auto y2 = X[lii + m] - X[rri + m];
				X[lri + m] = x1;
				X[lli + m] = x2;
				X[rri + m] = y1;
				X[rli + m] = y2;
			}

		}
	}
}

void fft_real(double* X, int N) {
	std::vector<int> N1s;
	std::vector<int> N2s;
	int M = N;
	while (M >= 4) {
		for (int n = 4; n <= M; n++) {
			const int n2 = n * n;
			while (M % n2 == 0) {
				N1s.push_back(n);
				M /= n2;
			}
			if (M % n == 0) {
				N2s.push_back(n);
				M /= n;
			}
		}
	}
	if (M % 2 == 0) {
		N2s.push_back(2);
	}
	if (M % 3 == 0) {
		N2s.push_back(3);
	}
	auto N3s = N1s;
	std::reverse(N3s.begin(), N3s.end());
	std::vector<int> Ns;
	Ns.insert(Ns.end(), N1s.begin(), N1s.end());
	Ns.insert(Ns.end(), N2s.begin(), N2s.end());
	Ns.insert(Ns.end(), N3s.begin(), N3s.end());
	const int N1 = std::reduce(N1s.begin(), N1s.end(), 1, std::multiplies<int>());
	const int N2 = std::reduce(N2s.begin(), N2s.end(), 1, std::multiplies<int>());
	const int N3 = N1;
	reverse_digits(X, Ns.data(), Ns.size(), 0, N1s.size());
	fft_real(X, N1, N1s.data(), N1s.data() + N1s.size(), N2 * N3);
	reverse_digits(X, Ns.data(), Ns.size(), N1s.size(), N1s.size() + N2s.size());
	fft_real(X, N2, N2s.data(), N2s.data() + N2s.size(), N3);
	if (N1 % 2 == 0) {
		auto rc = fft_pre_skew(X + N1 * N2 * N3 / 2, N2, N2s.data(), N2s.data() + N2s.size(), N3);
		fft_real(X + N1 * N2 * N3 / 2, N2, N2s.data(), N2s.data() + N2s.size(), N3);
		fft_post_skew(X + N1 * N2 * N3 / 2, N2, N3, rc);
	}
	for (int k1 = 1; k1 < N1; k1++) {
		fft_real(X + k1 * N2 * N3, N2, N2s.data(), N2s.data() + N2s.size(), N3);
	}
	transpose(X, N1, N2);
	reverse_digits(X, Ns.data(), Ns.size(), N1s.size() + N2s.size(), Ns.size());
}
