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

void reverse_digits(double* X, int* digs, int ndigs, int M) {
	const int N = std::reduce(digs, digs + ndigs, 1, std::multiplies<int>());
	for (int i = 0; i < N; i++) {
		int k = i;
		int j = 0;
		for (int l = 0; l < ndigs; l++) {
			std::div_t q = std::div(k, digs[l]);
			j *= digs[l];
			j += q.rem;
			k = q.quot;
		}
		if (i > j) {
			const int i0 = M * i;
			const int j0 = M * j;
			for (int m = 0; m < M; m++) {
				std::swap(X[i0 + m], X[j0 + m]);
			}
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

void reverse_imaginary(double* Y, int NHI, int NMID, int NLO) {
	for (int nhi = 0; nhi < NHI; nhi++) {
		for (int i = NMID; i < NMID - i - 1; i++) {
			const int k = NLO * (i + NMID * nhi);
			const int l = NLO * ((NMID - i - 1) + NMID * nhi);
			for (int nlo = 0; nlo < NLO; nlo++) {
				std::swap(Y[k + nlo], Y[l + nlo]);
			}
		}
	}
}

void fft_complex(double* X, double* Y, int N, int* digb, int* dige, int M) {
	const int N1 = *digb;
	const int N2 = N / N1;
	if (N2 > 1) {
		for (int n1 = 0; n1 < N1; n1++) {
			fft_complex(X + n1 * N2, Y + n1 * N2, N2, digb + 1, dige, M);
		}
	}
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
				sfft_complex(X + M * k2 + m, Y + M * k2 + m, M * N2, N1, wn.data(), w);
			}
		}
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

void fft_skew(double* X, int N, int* digb, int* dige, int M) {
	if (N % 2 == 0) {
	} else {

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
	reverse_digits(X, Ns.data(), Ns.size(), N1s.size(), N1s.size() + N2s.size());
	transpose(X, N1, N2);
	reverse_digits(X, Ns.data(), Ns.size(), N1s.size() + N2s.size(), Ns.size());
}
