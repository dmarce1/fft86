#include <fft86/fft.hpp>
#include <cstring>
#include <functional>
#include <vector>

void pfa(std::complex<double>* X, int N) {
	static thread_local std::vector<int> ifax;
	static thread_local std::vector<int> I;
	static thread_local std::vector<std::complex<double>> Y;
	ifax.resize(0);
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
		printf("%i ", nq);
		ifax.push_back(nq);
	}
	printf("\n");
	const int nfax = ifax.size();
	for (int k = 0; k < nfax; k++) {
		int R;
		int dN;
		const int N1 = ifax[k];
		const int N2 = N / N1;
		for (int j = 1; j <= N1; j++) {
			R = j;
			dN = j * N2;
			if (dN % N1 == 1) {
				break;
			}
		}
		Y.resize(N);
		for (int l = 0; l < N2; l++) {
			const int o = (l * (dN - 1)) % N;
			int k2 = o;
			for (int k = 0; k < N1; k++) {
				int n1 = o;
				Y[k2] = 0.0;
				for (int n = 0; n < N1; n++) {
					Y[k2] += X[n1] * std::polar(1.0, -2.0 * M_PI * n * k * R / N1);
					n1 = (n1 + dN) % N;
				}
				k2 = (k2 + dN) % N;
			}
		}
		for (int k = 0; k < N; k++) {
			X[k] = Y[k];
		}
	}
}
