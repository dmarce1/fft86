#include <fft86/fft.hpp>
#include <cstring>
#include <functional>
#include <vector>

void pfa(std::complex<double>* X, int N) {
	static thread_local std::vector<int> Ns;
	static thread_local std::vector<int> I;
	static thread_local std::vector<std::complex<double>> Y;
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
	for (int k = 0; k < Ns.size(); k++) {
		int R;
		int dN1, dN2;
		const int N1 = Ns[k];
		const int N2 = N / N1;
		for (int j = 1; j <= N1; j++) {
			R = j;
			dN1 = j * N2;
			if (dN1 % N1 == 1) {
				break;
			}
		}
		for (int j = 1; j <= N2; j++) {
			dN2 = j * N1;
			if (dN2 % N2 == 1) {
				break;
			}
		}
		Y.resize(N);
		for (int l = 0; l < N2; l++) {
			const int o = (l * dN2) % N;
			for (int k = 0; k < N1 / 2 + 1; k++) {
				const int k2r = (o + k * dN1) % N;
				const int k2i = (o + ((N1 - k) % N1) * dN1) % N;
				Y[k2r] = Y[k2i] = 0.0;
				for (int n = 0; n < N1; n++) {
					const int n1 = (o + n * dN1) % N;
					const auto w = std::polar(1.0, -2.0 * M_PI * n * k * R / N1);
					Y[k2r] += X[n1] * w.real();
					if (k2i != k2r) {
						Y[k2i] += X[n1] * w.imag();
					}
				}
			}
		}
		for (int k = 0; k < N; k++) {
			X[k] = Y[k];
		}
	}
	for (int k = Ns.size() - 1; k >= 0; k--) {
		int dN1, dN2;
		const int N1 = Ns[k];
		const int N2 = N / N1;
		for (int j = 1; j <= N1; j++) {
			dN1 = j * N2;
			if (dN1 % N1 == 1) {
				break;
			}
		}
		for (int j = 1; j <= N2; j++) {
			dN2 = j * N1;
			if (dN2 % N2 == 1) {
				break;
			}
		}
		Y.resize(N);
		std::complex<double> J(0, 1);
		for (int l = 0; l < N2; l++) {
			const int o = (l * dN2) % N;
			for (int k = 1; k < (N1 + 1) / 2; k++) {
				const int k2r = (o + k * dN1) % N;
				const int k2i = (o + ((N1 - k) % N1) * dN1) % N;
				const auto yr = Y[k2r] + J * Y[k2i];
				const auto yi = Y[k2r] - J * Y[k2i];
				Y[k2r] = yr;
				Y[k2i] = yi;
			}
			for (int k = 0; k < N; k++) {
				X[k] = Y[k];
			}
		}
	}
}
