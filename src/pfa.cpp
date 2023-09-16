#include <fft86/fft.hpp>
#include <cstring>
#include <functional>
#include <vector>

void pfa(double* X, int N) {
	static thread_local std::vector<int> Ns;
	static thread_local std::vector<int> I;
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
	const int nfax = Ns.size();
	for (int k = 0; k < nfax; k++) {
		int R;
		int dN, dM;
		const int N1 = Ns[k];
		const int N2 = N / N1;
		for (int j = 1; j <= N1; j++) {
			R = j;
			dN = j * N2;
			if (dN % N1 == 1) {
				break;
			}
		}
		for (int j = 1; j <= N1; j++) {
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
					const auto w = std::polar(1.0, -2.0 * M_PI * n * k * R / N1);
					Y[k2] += X[n1] * (w.real() - w.imag());
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
		printf("%e ", Wr[n]);
	}
	printf("\n");
	for (int n = 0; n < (1 << Ns.size()); n++) {
		printf("%e ", Wi[n]);
	}
	printf("\n");

}
