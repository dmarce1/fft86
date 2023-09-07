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

void fft_prime_factor_algorithm(double* X, const std::vector<int>& N) {
	std::vector<int> M(N.size());
	std::vector<int> L(N.size());
	int N0 = 1;
	for (int n = 0; n < N.size(); n++) {
		L[n] = N0;
		N0 *= N[n];
	}
	std::vector<double> Y(2 * N0);
	for (int n = 0; n < M.size(); n++) {
		M[n] = N0 / N[n];
	}
	std::vector<bool> visited(2 * N0, false);
	for (int n0 = 0; n0 < 2 * N0; n0++) {
		if (!visited[n0]) {
			int n = n0;
			int next;
			do {
				int p = n >> 1;
				for (int m = 0; m < N.size(); m++) {
					int d = p % N[m];
					p /= N[m];
					next += d * M[m];
				}
				next = next % N0;
				if (n & 1) {
					next += N0;
				}
				visited[next] = true;
				std::swap(X[n], X[next]);
				n = next;
			} while (next != n0);
		}
	}
	for (int i = 0; i < N.size(); i++) {
		int n1 = L[i];
		int o2 = L[i] * N[i];
		int n2 = N0 / L[i] / N[i];
		for (int j2 = 0; j2 < n2; j2++) {
			int j1 = 0;
			for (; j1 < n1 - 3; j1 += 4) {
				double* x = X + (o2 * j2 + j1);
				double* y = x + N0;
				sfft_complex_w4(x, y, L[i], N[i]);
			}
			for (; j1 < n1 - 1; j1 += 2) {
				double* x = X + (o2 * j2 + j1);
				double* y = x + N0;
				sfft_complex_w2(x, y, L[i], N[i]);
			}
			for (; j1 < n1; j1++) {
				double* x = X + (o2 * j2 + j1);
				double* y = x + N0;
				sfft_complex_w1(x, y, L[i], N[i]);
			}
		}
	}
	std::fill(visited.begin(), visited.end(), false);
	for (int n0 = 0; n0 < 2 * N0; n0++) {
		if (!visited[n0]) {
			int n = n0;
			int next;
			do {
				int next = 0;
				for (int m = 0; m < N.size(); m++) {
					next += L[m] * ((n >> 1) % N[m]);
				}
				visited[next] = true;
				std::swap(X[n], X[next]);
				n = next;
			} while (next != n0);
		}
	}

}

std::vector<std::pair<int, int>> compute_factors(int N) {
	int M = N;
	std::vector<int> p;
	std::vector<std::pair<int, int>> f;
	std::vector<std::pair<int, int>> factors;
	while (M > 1) {
		for (int m = 2; m <= M; m++) {
			if (M % m == 0) {
				int q = 0;
				do {
					q++;
					M /= m;
				} while (M % m == 0);
				if (q & 1) {
					p.push_back(m);
					q--;
				}
				if (q) {
					f.push_back(std::make_pair(m, q));
				}
			}
		}
	}
	for (auto fac : f) {
		factors.push_back(std::make_pair(fac.first, fac.second >> 1));
	}
	for (auto prm : p) {
		factors.push_back(std::make_pair(prm, 1));
	}
	std::reverse(f.begin(), f.end());
	for (auto fac : f) {
		factors.push_back(std::make_pair(fac.first, fac.second >> 1));
	}
	return std::move(factors);
}

void sort_mixed_radix(std::complex<double>* X, int N) {
	auto facs = compute_factors(N);
//	sort_mixed_radix(X, N, facs);

}

void sort_nonsquare(double* X, int* N, int nfacs, int NLO) {
	const int N0 = std::reduce(N, N + nfacs, 1, std::multiplies<int>());
	int next;
	int current = 1;
	do {
		int k = current;
		next = 0;
		for (int l = 0; l < nfacs; l++) {
			next *= N[l];
			next += k % N[l];
			k /= N[l];
		}
		for (int m = 0; m < NLO; m++) {
			std::swap(X[NLO * current + m], X[NLO * next + m]);
		}
		current = next;
	} while (next != 1);
}

void fft_sort(double* X, double* Y, int* Nhead, int* Ntail, int N, int NLO = 1) {
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
		//	printf( "%i %i\n", R, N1);
		const int NMID = N / (NLO * N1 * N1);
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
		N /= N1;
		NLO *= N1;
		for (int n1 = 0; n1 < N1; n1++) {
			fft_sort(X + n1 * N, Y + n1 * N, Nhead, Ntail, N, NLO);
		}
	} else {
		int* beg = Nhead;
		int* end = Ntail;
		static thread_local std::vector<bool> visited;
		static thread_local std::vector<double> tmpx;
		static thread_local std::vector<double> tmpy;
		tmpx.resize(NLO);
		tmpy.resize(NLO);
		N /= NLO;
		visited.resize(N);
		std::fill(visited.begin(), visited.end(), false);
		for (int n = 0; n < N; n++) {
			if (!visited[n]) {
				int current;
				int next = n;
				do {
					current = next;
					next = 0;
					int k = current;
					for (auto* n = end; n >= beg; n--) {
						next *= *n;
						next += k % *n;
						k /= *n;
					}
					if (current == n) {
						for (int m = 0; m < NLO; m++) {
							const int j = m + NLO * current;
							tmpx[m] = X[j];
							tmpy[m] = Y[j];
						}
					}
					for (int m = 0; m < NLO; m++) {
						const int i = m + NLO * next;
						std::swap(tmpx[m], X[i]);
						std::swap(tmpy[m], Y[i]);
					}
					visited[next] = true;
				} while (next != n);
			}
		}
	}
}

void fft_dif(double* X, double* Y, int N, int* Rptr, int M) {
	const auto* W0 = get_twiddles(N);
	static thread_local std::vector<std::complex<double>> w;
	const int N1 = *Rptr;
	const int N2 = N / N1;
	const int s = N2 * M;
	w.resize(N1);
	for (int n2 = 0; n2 < N2; n2++) {
		for (int n1 = 1; n1 < N1; n1++) {
			w[n1] = W0[n1 * n2];
		}
		for (int m = 0; m < M; m++) {
			sfft_complex_dif_w1(X + m + M * n2, Y + m + M * n2, s, N1, (double*) w.data());
		}
	}
	if (N2 > 1) {
		for (int k1 = 0; k1 < N1; k1++) {
			fft_dif(X + k1 * s, Y + k1 * s, N2, Rptr + 1, M);
		}
	}
}

void fft_dit(double* X, double* Y, int N, int* Rptr, int M) {
	const auto* W0 = get_twiddles(N);
	static thread_local std::vector<std::complex<double>> w;
	const int N1 = *Rptr;
	const int N2 = N / N1;
	const int s = N2 * M;
	w.resize(N1);
	if (N2 > 1) {
		for (int n1 = 0; n1 < N1; n1++) {
			fft_dit(X + n1 * s, Y + n1 * s, N2, Rptr + 1, M);
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		for (int n1 = 1; n1 < N1; n1++) {
			w[n1] = W0[n1 * k2];
		}
		for (int m = 0; m < M; m++) {
			sfft_complex_dit_w1(X + m + M * k2, Y + m + M * k2, s, N1, (double*) w.data());
		}
	}
}

void apply_twiddles(double* X, double* Y, const std::vector<int>& N2s, int N1) {
	const int N2 = std::reduce(N2s.begin(), N2s.end(), 1, std::multiplies<int>());
	const int N = N1 * N2;
	const auto* W = get_twiddles(N);
	for (int k2 = 0; k2 < N2; k2++) {
		int k = k2;
		int k2r = 0;
		for (int n = 0; n < N2s.size(); n++) {
			k2r *= N2s[n];
			k2r += k % N2s[n];
			k /= N2s[n];
		}
		for (int n1 = 0; n1 < N1; n1++) {
			const auto w = W[n1 * k2r];
			auto& x = X[k2 * N1 + n1];
			auto& y = Y[k2 * N1 + n1];
			auto tmp = x;
			x = tmp * w.real() - y * w.imag();
			y = tmp * w.imag() + y * w.real();
		}
	}

}

void fft_complex(double* X, double* Y, int N) {
	std::vector<int> evens;
	std::vector<int> odds;
	std::vector<int> N1s;
	std::vector<int> N2s;
	int M = N;
	while (M > 1) {
		for (int n = 2; n <= M; n++) {
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
	fft_dif(X, Y, N1, N1s.data(), N2);
	std::reverse(N1s.begin(), N1s.end());
	apply_twiddles(X, Y, N1s, N2);
	fft_sort(X, Y, Ns.data(), Ns.data() + Ns.size() - 1, N1 * N2);
	std::reverse(N2s.begin(), N2s.end());
	fft_dit(X, Y, N2, N2s.data(), N1);

}

