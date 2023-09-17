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
	dNs.resize(Ns.size());
	for (int n = 0; n < Ns.size(); n++) {
		for (int j = 1; j <= N / Ns[n]; j++) {
			dNs[n] = j * (N / Ns[n]);
			if (dNs[n] % Ns[n] == 1) {
				break;
			}
		}
	}
	for (int KK = 0; KK < Ns.size(); KK++) {
		int R;
		int dN, dM;
		const int N1 = Ns[KK];
		const auto* W = get_twiddles(N1);
		const int N2 = N / N1;
		dN = dNs[KK];
		R = dN / N2;
		for (int j = 1; j <= N2; j++) {
			dM = j * N1;
			if (dM % N2 == 1) {
				break;
			}
		}
		Y.resize(N);
		int beg = 0;
		for (int l = 0; l < N2; l++) {
			int k2 = beg;
			for (int k = 0; k < N1; k++) {
				int n1 = beg;
				Y[k2] = 0.0;
				for (int n = 0; n < N1; n++) {
					const int Rnk = (R * n * k) % N1;
					Y[k2] += X[n1] * (W[Rnk].real() - W[Rnk].imag());
					n1 += dN;
					while (n1 >= N) {
						n1 -= N;
					}
				}
				k2 += dN;
				while (k2 >= N) {
					k2 -= N;
				}
			}
			beg += dM;
			while (beg >= N) {
				beg -= N;
			}
		}
		for (int k = 0; k < N; k++) {
			X[k] = Y[k];
		}
	}
	std::vector<double> Wr(1 << Ns.size(), 0);
	std::vector<double> Wi(1 << Ns.size(), 0);
	std::vector<double> qr(1 << Ns.size(), 0);
	std::vector<double> qi(1 << Ns.size(), 0);
	Wr[0] = 1;
	Wr[1] = 1;
	Wi[0] = -1;
	Wi[1] = 1;
	for (int n = 1; n < Ns.size(); n++) {
		int M = 1 << n;
		qr = Wr;
		qi = Wi;
		std::fill(Wr.begin(), Wr.end(), 0);
		std::fill(Wi.begin(), Wi.end(), 0);
		for (int m = 0; m < M; m++) {
			Wr[m] += qr[m] + qi[m];
			Wr[m + M] += qr[m] - qi[m];
			Wi[m] -= qr[m] - qi[m];
			Wi[m + M] += qi[m] + qr[m];
		}
	}
	double wtot = 0.0;
	for (int n = 0; n < (1 << Ns.size()); n++) {
		wtot += Wr[n];
	}
	wtot = 1.0 / wtot;
	for (int n = 0; n < (1 << Ns.size()); n++) {
		Wr[n] *= wtot;
		Wi[n] *= wtot;
	}
	std::vector<int> B(Ns.size());
	B[0] = 1;
	for (int n = 1; n < Ns.size(); n++) {
		B[n] = B[n - 1] * Ns[n - 1];
	}
	std::vector<int> J(1 << Ns.size(), 0);
	std::vector<int> digs(Ns.size(), 0);
	for (int k = 0; k < N; k++) {
		X[k] = Y[k];
	}
	for (int n = 0; n < N; n++) {
		int m = n;
		int q = 0;
		int nn = 0;
		for (int k = Ns.size() - 1; k >= 0; k--) {
			digs[k] = m / B[k];
			m = m % B[k];
			nn += digs[k] * dNs[k];
			nn = nn % N;
		}
		for (int j = 0; j < J.size(); j++) {
			auto D = digs;
			int m = j;
			int t = Ns.size() - 1;
			while (m) {
				if (m & 1) {
					D[t] = (Ns[t] - D[t]) % Ns[t];
				}
				m >>= 1;
				t--;
			}
			int k = 0;
			for (int q = Ns.size() - 1; q >= 0; q--) {
				k += D[q] * dNs[q];
				k = k % N;
			}
			J[j] = k;
		}
		Y[nn] = 0.0;
		for (int j = 0; j < J.size(); j++) {
			if (nn >= N / 2 + 1) {
				Y[nn] -= Wi[j] * X[J[j]];
			} else {
				Y[nn] += Wr[j] * X[J[j]];
			}
		}
	}
	for (int k = 0; k < N; k++) {
		X[k] = Y[k];
	}
}
