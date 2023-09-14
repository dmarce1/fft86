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

void sfft_real_complex_dit(double* X, double* Y, int s, int N1, double* w, int width) {
	switch (width) {
	case 4:
		return sfft_real_complex_dit_w4(X, Y, s, N1, w);
	case 2:
		return sfft_real_complex_dit_w2(X, Y, s, N1, w);
	case 1:
		return sfft_real_complex_dit_w1(X, Y, s, N1, w);
	}
}

void sfft_real_complex_dif(double* X, double* Y, int s, int N1, double* w, int width) {
	switch (width) {
	case 4:
		return sfft_real_complex_dif_w4(X, Y, s, N1, w);
	case 2:
		return sfft_real_complex_dif_w2(X, Y, s, N1, w);
	case 1:
		return sfft_real_complex_dif_w1(X, Y, s, N1, w);
	}
}

void fft_real(double* X, int k2, int* dbegin, int* dhead, int* dtail, int* dend) {
	if (dhead >= dtail) {
		return;
	}
	const int Nhi = std::reduce(dbegin, dhead, 1, std::multiplies<int>());
	const int Nmid = std::reduce(dhead, dtail, 1, std::multiplies<int>());
	const int Nlo = std::reduce(dtail, dend, 1, std::multiplies<int>());
	const int N = Nhi * Nmid * Nlo;
	int N1 = *dhead;
	int N2 = N / N1;
	if (k2 == 0) {
		for (int w = 4, n = 0; w >= 1; w >>= 1) {
			for (; n <= N2 - w; n++) {
				sfft_real(X + n, N2, N1, w);
			}
		}
	} else if (k2 == Nhi / 2 && (Nhi % 2 == 0)) {
		for (int w = 4, n = 0; w >= 1; w >>= 1) {
			for (; n <= N2 - w; n++) {
				sfft_skew(X + n, N2, N1, w);
			}
		}
	} else if (k2 < Nhi) {
		const auto* W0 = get_twiddles(N);
		std::array<std::vector<double>, 5> tw;
		for (int w = 1; w <= 4; w *= 2) {
			tw[w].resize(N1);
			for (int n1 = 1; n1 < N1; n1++) {
				for (int l = 0; l < w; l++) {
					tw[w][2 * n1 + l] = W0[n1 * k2].real();
					tw[w][2 * n1 + w + l] = W0[n1 * k2].imag();
				}
			}
		}
		for (int w = 4, n = 0; w >= 1; w /= 2) {
			for (; n <= N2 - w; n++) {
		//		sfft_complex(X + n, X + N + n, N2, N1, tw[w].data(), w);
			}
		}
	}
	for (int n1 = 0; n1 < N1; n1++) {
		fft_real(X + n1 * N2, Nhi * n1 + k2, dbegin, dhead + 1, dtail - 1, dend);
	}
//	for (int k1 = 1; k1 < N - k1; k1 += 2) {
//		std::swap(X[k1], X[N - k1]);
//	}
	if (*dhead == *dtail) {
		N1 = *dhead;
		const int Mhi = Nmid / *dhead;
		const int Mlo = Nlo / *dtail;
		for (int nhi = 0; nhi < Mhi; nhi++) {
			for (int n1 = 0; n1 < N1; n1++) {
				for (int n2 = 0; n2 < N1; n2++) {
					const int i = Mlo * (n2 + N1 * (nhi + Mhi * n1));
					const int j = Mlo * (n1 + N1 * (nhi + Mhi * n2));
					if (i < j) {
						for (int nlo = 0; nlo < Mlo; nlo++) {
							std::swap(X[i + nlo], X[j + nlo]);
						}
					}
				}
			}
		}
	} else {
		static thread_local std::vector<bool> visited;
		const int N1 = *dhead;
		const int N2 = Nmid / *dhead;
		const int N3 = *dtail;
		const int M = Nlo / *dtail;
		visited.resize(N1 * N3);
		std::fill(visited.begin(), visited.end(), false);
		for (int i = 1; i < N1 * N3 - 1; i++) {
			if (!visited[i]) {
				int j = i;
				do {
					j = N2 * N3 * (j % N1) + N3 * ((j / N1) % N2) + j / (N1 * N2);
					for (int m = 0; m < M; m++) {
						std::swap(X[m + M * j], X[m + M * i]);
					}
					visited[j] = true;
				} while (j != i);
			}
		}
	}
	const int& k2hi = k2;
	if (Nhi != 1) {
		N1 = *dtail;
		N2 = N / N1;
		for (int k2lo = 0; k2lo < N2 / 2 + 1; k2lo++) {
			const int k = k2hi * N2 + k2lo;
			if (k == 0) {
				for (int w = 4, n = 0; w >= 1; w >>= 1) {
					for (; n <= Nlo - w; n++) {
						sfft_real(X + Nlo * k2lo + n, N2, N1, w);
					}
				}
			} else if (k == N2 / 2 && (N2 % 2 == 0)) {
				for (int w = 4, n = 0; w >= 1; w >>= 1) {
					for (; n <= Nlo - w; n++) {
						sfft_skew(X + Nlo * k2lo + n, N2, N1, w);
					}
				}
			} else {
				const auto* W0 = get_twiddles(N);
				std::array<std::vector<double>, 5> tw;
				for (int w = 1; w <= 4; w *= 2) {
					tw[w].resize(N1);
					for (int n1 = 1; n1 < N1; n1++) {
						for (int l = 0; l < w; l++) {
							tw[w][2 * n1 + l] = W0[n1 * k2].real();
							tw[w][2 * n1 + w + l] = W0[n1 * k2].imag();
						}
					}
				}
				for (int w = 4, n = 0; w >= 1; w /= 2) {
					for (; n <= Nlo - w; n++) {
			//			sfft_complex(X + Nlo * k2lo + n, X + N + Nlo * k2lo + n, N2, N1, tw[w].data(), w);
					}
				}
			}
		}
	} else {
		sfft_real(X, N2, N1, 1);
		if (N2 % 2 == 0) {
			sfft_skew(X + N2 / 2, N2, N1, 1);
		}
		const auto* W0 = get_twiddles(N);
		std::array<std::vector<double>, 5> tw;
		for (int w = 4, n = 0; w >= 1; w /= 2) {
			for (int k2 = 1; k2 <= N2 / 2 - w; k2++) {
				for (int w = 1; w <= 4; w *= 2) {
					tw[w].resize(N1);
					for (int n1 = 1; n1 < N1; n1++) {
						for (int dk = 0; dk < w; dk++) {
							tw[w][2 * n1 + dk] = W0[n1 * (k2 + dk)].real();
							tw[w][2 * n1 + w + dk] = W0[n1 * (k2 + dk)].imag();
						}
					}
				}
			//	sfft_complex(X + k2, X + N + k2, N2, N1, tw[w].data(), w);
			}
		}
	}
}

void fft_real(double* X, int N) {
	std::vector<int> N1s;
	std::vector<int> N2s;
	std::vector<int> N3s;
	std::vector<int> Ns;
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
	N3s = N1s;
	std::reverse(N3s.begin(), N3s.end());
	Ns.insert(Ns.end(), N1s.begin(), N1s.end());
	Ns.insert(Ns.end(), N2s.begin(), N2s.end());
	Ns.insert(Ns.end(), N3s.begin(), N3s.end());
	fft_real(X, 0, Ns.data(), Ns.data(), Ns.data() + Ns.size(), Ns.data() + Ns.size());
}
