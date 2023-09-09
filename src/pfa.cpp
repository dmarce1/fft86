#include <fft86/fft.hpp>
#include <cstring>
#include <functional>
#include <vector>

void fft(std::complex<double>* X, int N, std::function<int(int)> P) {
	std::vector<std::complex<double>> Y(N);
	for (int k = 0; k < N; k++) {
		Y[P(k)] = 0.0;
		for (int n = 0; n < N; n++) {
			const auto W = std::polar(1.0, -2.0 * M_PI * n * k / N);
			Y[P(k)] += X[P(n)] * W;
		}
	}
	std::memcpy(X, Y.data(), sizeof(std::complex<double>) * N);
}

void fft_pfa(std::complex<double>* X, int N1, int N2) {
	const int N = N1 * N2;
	std::vector<std::complex<double>> Y(N);
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			Y[n1 * N2 + n2] = X[(n1 * N2 + n2 * N1) % N];
		}
	}
	for (int n1 = 0; n1 < N1; n1++) {
		fft(Y.data() + n1 * N2, N2, std::function<int(int)>([N1,N2,n1](int n2) {
			return (n1 * N2 + n2 * N1) % (N1*N2);
		}));
	}
	for (int n1 = 0; n1 < N1; n1++) {
		for (int k2 = 0; k2 < N2; k2++) {
			X[k2 * N1 + n1] = Y[n1 * N2 + k2];
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		fft(X + k2 * N1, N1, std::function<int(int)>([N1,N2,k2](int n1) {
			return (n1 * N1 + k2 * N2) % (N1*N2);
		}));
	}
	for (int k = 0; k < N; k++) {
		const int k1 = k % N1;
		const int k2 = k % N2;
		Y[k] = X[k2 * N1 + k1];
	}
	std::memcpy(X, Y.data(), sizeof(std::complex<double>) * N);
}
