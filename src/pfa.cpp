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
		int mu;
		int MM;
		const int ifac = ifax[k];
		const int M = N / ifac;
		for (int j = 1; j <= ifac; j++) {
			mu = j;
			MM = j * M;
			if (MM % ifac == 1) {
				break;
			}
		}
		printf("ifac: %i M: %i MM: %i mu: %i\n", ifac, M, MM, mu);
		I.resize(ifac);
		Y.resize(N);
		I[0] = 0;
		for (int m = 0; m < ifac - 1; m++) {
			I[m + 1] = (I[m] + MM) % N;
		}
		for (int l = 0; l < M; l++) {
			for (int m = 0; m < ifac; m++) {
				printf("%i\n", I[m]);
			}
			printf("\n");
			for (int k = 0; k < ifac; k++) {
				Y[I[k]] = 0.0;
				for (int n = 0; n < ifac; n++) {
					Y[I[k]] += X[I[n]] * std::polar(1.0, -2.0 * M_PI * n * k * mu / ifac);
				}
			}
			const int ix = I[ifac - 1] + 1;
			for (int m = ifac - 2; m >= 0; m--) {
				I[m + 1] = I[m] + 1;
			}
			I[0] = ix;
		}
		printf("\n");
		for (int k = 0; k < N; k++) {
			X[k] = Y[k];
		}
	}
}
