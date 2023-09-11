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


void fft_sort(double* X, double* Y, int* Nhead, int* Ntail, int N, int NLO = 1) {
	static thread_local std::vector<bool> visited;
	if (Nhead >= Ntail) {
		return;
	}
	int N1;
	int R;
	if (*Nhead == *Ntail) {
		N1 = 1;
		R = *Nhead;
		while (*Nhead == R && *Ntail == R) {
			N1 *= R;
			Nhead++;
			Ntail--;
			if (Nhead >= Ntail) {
				break;
			}
		}
		const int NMID = N / (NLO * N1 * N1);
		if (NLO == 1) {
			scramble_hi(X, R, N1, N / N1);
			scramble_hi(Y, R, N1, N / N1);
			transpose(X, N1, N / (N1 * N1));
			transpose(Y, N1, N / (N1 * N1));
			scramble_hi(X, R, N1, N / N1);
			scramble_hi(Y, R, N1, N / N1);
		} else {
			const int Rm1 = R - 1;
			int n1r = 0;
			for (int n1 = 0; n1 < N1; n1++) {
				for (int n2 = 0; n2 < NMID; n2++) {
					int n3r = 0;
					for (int n3 = 0; n3 < N1; n3++) {
						int i = n3 + N1 * (n2 + NMID * n1);
						int j = n1r + N1 * (n2 + NMID * n3r);
						if (i < j) {
							for (int n4 = 0; n4 < NLO; n4++) {
								const int i0 = i * NLO + n4;
								const int j0 = j * NLO + n4;
								std::swap(X[i0], X[j0]);
								std::swap(Y[i0], Y[j0]);
							}
						}
						if (n3r != N1 - 1) {
							int k = N1 / R;
							n3r++;
							while (Rm1 * k < n3r) {
								n3r -= Rm1 * k;
								k /= R;
							}
							n3r += k - 1;
						}
					}
				}
				if (n1r != N1 - 1) {
					int k = N1 / R;
					n1r++;
					while (Rm1 * k < n1r) {
						n1r -= Rm1 * k;
						k /= R;
					}
					n1r += k - 1;
				}
			}
		}
		N /= N1;
		NLO *= N1;
		for (int n1 = 0; n1 < N1; n1++) {
			fft_sort(X + n1 * N, Y + n1 * N, Nhead, Ntail, N, NLO);
		}
	} else {
		int* beg = Nhead;
		int* end = Ntail;
		while (beg != end) {
			int N0 = N / NLO;
			visited.resize(N0);
			std::fill(visited.begin(), visited.end(), false);
			int N2 = 1;
			int N1 = *beg;
			for (auto i = beg + 1; i <= end; i++) {
				N2 *= *i;
			}
			const int NNm1 = N1 * N2 - 1;
			for (int n = 1; n < NNm1; n++) {
				if (!visited[n]) {
					int current;
					int next = n;
					const int j0 = NLO * n;
					do {
						current = next;
						next = (current * N1) % NNm1;
						visited[next] = true;
						if (next != n) {
							const int i0 = NLO * next;
							for (int m = 0; m < NLO; m++) {
								const int i = m + i0;
								const int j = m + j0;
								std::swap(X[j], X[i]);
								std::swap(Y[j], Y[i]);
							}
						}
					} while (next != n);
				}
			}
			beg++;
			NLO *= N1;
		}
	}
}

void fft_complex_dif(double* X, double* Y, int N, int* Rptr, int M) {
	const auto* W0 = get_twiddles(N);
	static thread_local std::vector<double> w;
	const int N1 = *Rptr;
	const int N2 = N / N1;
	const int s = N2 * M;
	w.resize(8 * N1);
	for (int n2 = 0; n2 < N2; n2++) {
		int m = 0;
		for (int n1 = 1; n1 < N1; n1++) {
			for (int l = 0; l < 4; l++) {
				w[8 * n1 + l] = W0[n1 * n2].real();
				w[8 * n1 + l + 4] = W0[n1 * n2].imag();
			}
		}
		for (; m < M - 3; m += 4) {
			sfft_complex_dif_w4(X + m + M * n2, Y + m + M * n2, s, N1, w.data());
		}
		for (int n1 = 1; n1 < N1; n1++) {
			for (int l = 0; l < 2; l++) {
				w[4 * n1 + l] = W0[n1 * n2].real();
				w[4 * n1 + l + 2] = W0[n1 * n2].imag();
			}
		}
		for (; m < M - 1; m += 2) {
			sfft_complex_dif_w2(X + m + M * n2, Y + m + M * n2, s, N1, w.data());
		}
		for (int n1 = 1; n1 < N1; n1++) {
			w[2 * n1] = W0[n1 * n2].real();
			w[2 * n1 + 1] = W0[n1 * n2].imag();
		}
		for (; m < M; m++) {
			sfft_complex_dif_w1(X + m + M * n2, Y + m + M * n2, s, N1, w.data());
		}
	}
	if (N2 > 1) {
		for (int k1 = 0; k1 < N1; k1++) {
			fft_complex_dif(X + k1 * s, Y + k1 * s, N2, Rptr + 1, M);
		}
	}
}

void fft_real_dif(double* X, int N, int* Rptr, int M);

void fft_skew_dif(double* X, int N, int* Rptr, int M) {
	/*const auto* W2 = get_twiddles(2 * N);
	if (N % 2 == 0) {
		bool is2;
		const double xny = X[N / 2];
		X[N / 2] = 0.0;
		double* xp = X;
		double* xm = X + N / 2;
		for (int n = 1; n < N; n++) {
			auto tmp = X[n];
			X[n] = tmp - X[N - n];
			X[N - n] = -tmp - X[N - n];
		}
		if ((N / 2) % 2 == 0) {
			xp[N / 4] *= sqrt(0.5);
			xm[N / 4] *= sqrt(0.5);
		}
		for (int n = 1; n < N / 2 - n; n++) {
			const double c = W2[2 * N - n].real();
			const double s = W2[2 * N - n].imag();
			auto tmpp = xp[n];
			xp[n] = 0.5 * (tmpp * (c + s) + xp[N / 2 - n] * (s - c));
			xp[N / 2 - n] = 0.5 * (tmpp * (c - s) + xp[N / 2 - n] * (s + c));
			auto tmpm = xm[n];
			xm[n] = 0.5 * (tmpm * (c + s) + xm[N / 2 - n] * (s - c));
			xm[N / 2 - n] = 0.5 * (tmpm * (c - s) + xm[N / 2 - n] * (s + c));
		}
		const bool is2 = (*Rptr == 2);
		if (is2) {
			Rptr++;
		} else {
			*Rptr >>= 1;
		}
		fft_real_dif(xp, N / 2, Rptr, M);
		if (!is2) {
			*Rptr <<= 1;
		}
		const int d = 1 - 2 * (N % 2);
		for (int n = 1; n < (N + 2) / 4; n++) {
			auto tmpp = xp[2 * n];
			auto tmpm = xm[2 * n];
			xp[2 * n] = tmpp + xp[2 * n + d];
			xm[2 * n] = tmpm + xm[2 * n + d];
			xp[2 * n + d] = tmpp - xp[2 * n + d];
			xm[2 * n + d] = tmpm - xm[2 * n + d];
		}
		for (int n = 1; n < (N + 2) / 4; n++) {
			X[2 * n + d] = std::pow(-1, n / 2) * (X[2 * n + d] - xny);
		}
	} else {
		const int dn = N / 2;
		for (int n = 1; n < N; n += 2) {
			X[n] = -X[n];
		}
		fft_real_dif(X, N, Rptr, M);
		std::complex<double> z0(0.0, 0.0);
		for (int n = 0; n < N; n++) {
			//X[n] *= W2[n];
		}
		int n = 0;
		do {
			int l = n;
			n += (N - 1) / 2;
			n %= N;
			std::swap(X[n], X[l]);
		} while (n != 0);
		for (int n = 0; n < N / 2; n++) {
			X[2 * n + 1] = -X[2 * n + 1];
		}
		for (int n = 1; n < N / 2; n++) {
			std::swap(X[N - n], X[n]);
		}
		X[0] = z0.real();
		X[1] = z0.imag();
	}*/
}

void fft_real_dif(double* X, int N, int* Rptr, int M) {
	const int N1 = *Rptr;
	const int N2 = N / N1;
	const int s = N2 * M;
	for (int n2 = 0; n2 < N2; n2++) {
		int m = 0;
		for (; m < M - 3; m += 4) {
			sfft_real_w4(X + m + M * n2, s, N1);
		}
		for (; m < M - 1; m += 2) {
			sfft_real_w2(X + m + M * n2, s, N1);
		}
		for (; m < M; m++) {
			sfft_real_w1(X + m + M * n2, s, N1);
		}
	}
	if (N2 > 1) {
		fft_real_dif(X, N2, Rptr + 1, M);
		if (N1 % 2 == 0) {
			fft_skew_dif(X, N2, Rptr + 1, M);
		}
		const int ks = 2 - (N1 % 2);
		const int dk = 1 - (N1 % 2);
		for (int k1 = ks; k1 < N1; k1++) {
			auto* x = X + k1 * s;
			auto* y = X + (N1 - k1 + dk) * s;
			fft_complex_dif(x, y, N2, Rptr + 1, M);
		}
	}
}

void fft_complex_dit(double* X, double* Y, int N, int* Rptr, int M) {
	const auto* W0 = get_twiddles(N);
	static thread_local std::vector<double> w;
	const int N1 = *Rptr;
	const int N2 = N / N1;
	const int s = N2 * M;
	w.resize(8 * N1);
	if (N2 > 1) {
		for (int n1 = 0; n1 < N1; n1++) {
			fft_complex_dit(X + n1 * s, Y + n1 * s, N2, Rptr + 1, M);
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		int m = 0;
		for (int n1 = 1; n1 < N1; n1++) {
			for (int l = 0; l < 4; l++) {
				w[8 * n1 + l] = W0[n1 * k2].real();
				w[8 * n1 + l + 4] = W0[n1 * k2].imag();
			}
		}
		for (; m < M - 3; m += 4) {
			sfft_complex_dit_w4(X + m + M * k2, Y + m + M * k2, s, N1, w.data());
		}
		for (int n1 = 1; n1 < N1; n1++) {
			for (int l = 0; l < 2; l++) {
				w[4 * n1 + l] = W0[n1 * k2].real();
				w[4 * n1 + l + 2] = W0[n1 * k2].imag();
			}
		}
		for (; m < M - 1; m += 2) {
			sfft_complex_dit_w2(X + m + M * k2, Y + m + M * k2, s, N1, w.data());
		}
		for (int n1 = 1; n1 < N1; n1++) {
			w[2 * n1] = W0[n1 * k2].real();
			w[2 * n1 + 1] = W0[n1 * k2].imag();
		}
		for (; m < M; m++) {
			sfft_complex_dit_w1(X + m + M * k2, Y + m + M * k2, s, N1, w.data());
		}
	}
}

void apply_twiddles(double* X, double* Y, const std::vector<int>& N2s, int N1) {
	const int N2 = std::reduce(N2s.begin(), N2s.end(), 1, std::multiplies<int>());
	const int N = N1 * N2;
	const auto& C = *get_6step_cos_twiddles(N1, N2);
	const auto& S = *get_6step_sin_twiddles(N1, N2);
	const auto* W = get_twiddles(N);
	for (int k2 = 0; k2 < N2; k2++) {
		int k = k2;
		int k2r = 0;
		for (int n = 0; n < N2s.size(); n++) {
			k2r *= N2s[n];
			k2r += k % N2s[n];
			k /= N2s[n];
		}
		int n1 = 0;
		for (; n1 < N1 - 3; n1 += 4) {
			v4df x, y, c, s;
			const int i = k2 * N1 + n1;
			const int j = k2r * N1 + n1;
			x.load(X + i);
			y.load(Y + i);
			c.load(C.data() + j);
			s.load(S.data() + j);
			auto tmp = x;
			x = tmp * c - y * s;
			y = tmp * s + y * c;
			x.store(X + i);
			y.store(Y + i);
		}
		for (; n1 < N1 - 1; n1 += 2) {
			v2df x, y, c, s;
			const int i = k2 * N1 + n1;
			const int j = k2r * N1 + n1;
			x.load(X + i);
			y.load(Y + i);
			c.load(C.data() + j);
			s.load(S.data() + j);
			auto tmp = x;
			x = tmp * c - y * s;
			y = tmp * s + y * c;
			x.store(X + i);
			y.store(Y + i);
		}
		for (; n1 < N1; n1++) {
			v1df x, y, c, s;
			const int i = k2 * N1 + n1;
			const int j = k2r * N1 + n1;
			x.load(X + i);
			y.load(Y + i);
			c.load(C.data() + j);
			s.load(S.data() + j);
			auto tmp = x;
			x = tmp * c - y * s;
			y = tmp * s + y * c;
			x.store(X + i);
			y.store(Y + i);
		}
	}

}

static timer tm1, tm2, tm3, tm4, tm5;
void fft_complex(double* X, double* Y, int N) {

	tm1.start();
	std::vector<int> evens;
	std::vector<int> odds;
	std::vector<int> N1s;
	std::vector<int> N2s;
	int M = N;
	while (M >= 4) {
		for (int n = 4; n <= M; n++) {
			const int n2 = n * n;
			while (M % n2 == 0) {
				evens.push_back(n);
				M /= n2;
			}
			if (M % n == 0) {
				odds.push_back(n);
				M /= n;
			}
		}
	}
	if (M % 2 == 0) {
		odds.push_back(2);
	}
	if (M % 3 == 0) {
		odds.push_back(3);
	}
//	std::reverse(evens.begin(), evens.end());
	//std::reverse(odds.begin(), odds.end());
	N1s = evens;
	M = std::reduce(odds.begin(), odds.end(), 1, std::multiplies<int>());
	int i = 0;
	while (M > 1) {
		N2s.push_back(odds[i]);
		M /= odds[i];
		i++;
		if (M > 1) {
			N1s.push_back(odds[i]);
			M /= odds[i];
			i++;
		}
	}
	std::reverse(evens.begin(), evens.end());
	N2s.insert(N2s.end(), evens.begin(), evens.end());
	const int N2 = std::reduce(N2s.begin(), N2s.end(), 1, std::multiplies<int>());
	const int N1 = std::reduce(N1s.begin(), N1s.end(), 1, std::multiplies<int>());
	std::vector<int> Ns = N1s;
	Ns.insert(Ns.end(), N2s.begin(), N2s.end());
	/*	printf("N1 = ");
	 for (int n = 0; n < N1s.size(); n++) {
	 printf("(%i)", N1s[n]);
	 }
	 printf("\n");
	 printf("N2 = ");
	 for (int n = 0; n < N2s.size(); n++) {
	 printf("(%i)", N2s[n]);
	 }
	 printf("\n");*/
	tm1.stop();
	tm2.start();
	fft_complex_dif(X, Y, N1, N1s.data(), N2);
	tm2.stop();
	std::reverse(N1s.begin(), N1s.end());
	tm3.start();
	apply_twiddles(X, Y, N1s, N2);
	tm3.stop();
	tm4.start();
	fft_sort(X, Y, Ns.data(), Ns.data() + Ns.size() - 1, N1 * N2);
	tm4.stop();
	std::reverse(N2s.begin(), N2s.end());
	tm5.start();
	fft_complex_dit(X, Y, N2, N2s.data(), N1);
	tm5.stop();
	const double tinv = 100.0 / (tm1.read() + tm2.read() + tm3.read() + tm4.read() + tm5.read());
//	printf("%e %e %e %e %e\n", tm1.read() * tinv, tm2.read() * tinv, tm3.read() * tinv, tm4.read() * tinv,
//			tm5.read() * tinv);

}

