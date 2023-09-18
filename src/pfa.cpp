#include <fft86/fft.hpp>
#include <fft86/twiddles.hpp>
#include <cstring>
#include <functional>
#include <vector>

void pfa(double* X, int N) {
	static thread_local std::vector<int> Ns;
	static thread_local std::vector<int> dNs;
	static thread_local std::vector<double> Y;
	Ns.resize(0);
	int M = N;
	int n1 = 2;
	while (M > 1) {
		int nq;
		for (int n = n1; n <= M; n++) {
			if (M % n == 0) {
				n1 = n;
				nq = 1;
				while (M % n == 0) {
					nq *= n;
					M /= n;
				}
				break;
			}
		}
		Ns.push_back(nq);
	}
	int Nhi, Nlo;
	Nhi = N;
	Nlo = 1;
	int dNhi = 0;
	int dNlo = 0;
	for (int KK = 0; KK < Ns.size(); KK++) {
		int R;
		int dN1;
		const int N1 = Ns[KK];
		const auto* W = get_twiddles(N1);
		Nhi /= N1;
		for (int j = 1;; j++) {
			dN1 = j * (N / N1);
			R = j;
			if (dN1 % N1 == 1) {
				break;
			}
		}
		if (Nlo > 1) {
			for (int j = 1;; j++) {
				dNlo = j * (N / Nlo);
				if (dNlo % Nlo == 1) {
					break;
				}
			}
		}
		if (Nhi > 1) {
			for (int j = 1;; j++) {
				dNhi = j * (N / Nhi);
				if (dNhi % Nhi == 1) {
					break;
				}
			}
		}
		Y.resize(N);
		for (int nhi = 0; nhi < Nhi; nhi++) {
			int k2hi, k2lo;
			k2hi = (nhi * dNhi) % N;
			k2lo = k2hi;
			for (int nlo = 0; nlo < Nlo; nlo++) {
				int k2lo0 = k2lo;
				int k2hi0 = k2hi;
				for (int k = 0; k < N1; k++) {
					const int k2lo = (nhi * dNhi + nlo * dNlo + k * dN1) % N;
					const int k2hi = (nhi * dNhi + nlo * dNlo + ((N1 - k) % N1) * dN1) % N;
					if (k < N1 / 2 + 1) {
						Y[k2hi] = 0.0;
						Y[k2lo] = 0.0;
						for (int n = 0; n < N1; n++) {
							const int n1lo = (nhi * dNhi + nlo * dNlo + n * dN1) % N;
							const int n1hi = (nhi * dNhi + ((Nlo - nlo) % Nlo) * dNlo + ((N1 - n) % N1) * dN1) % N;
							const int Rnk = (R * n * k) % N1;
							if (KK == 0) {
								Y[k2lo] += X[n1lo] * (W[Rnk].real() - W[Rnk].imag());
								if (k2hi != k2lo) {
									Y[k2hi] += X[n1lo] * (W[Rnk].real() + W[Rnk].imag());
								}
							} else {
								Y[k2lo] += X[n1lo] * W[Rnk].real() + X[n1hi] * W[Rnk].imag();
								if (k2hi != k2lo) {
									Y[k2hi] += X[n1lo] * W[Rnk].real() - X[n1hi] * W[Rnk].imag();
								}
							}
						}
					}
				}
			}
		}
		for (int k = 0; k < N; k++) {
			X[k] = Y[k];
		}
		Nlo *= N1;
	}
	for (int n = 0; n <= N - n; n++) {
		const int p = (N - n) % N;
		const auto hp = X[p];
		const auto hn = X[n];
		X[p] = 0.5 * (hp - hn);
		X[n] = 0.5 * (hp + hn);
	}
}
