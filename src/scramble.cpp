#include <fft86/defs.hpp>
#include <fft86/scramble.hpp>

#include <algorithm>
#include <cassert>
#include <cstring>
#include <numeric>
#include <map>
#include <unordered_map>

void scramble_hi(double* X, int R, int N1, int NLO) {
	int nr = 0;
	int Rm1 = R - 1;
	for (int n1 = 0; n1 < N1; n1++) {
		if (n1 < nr) {
			for (int n2 = 0; n2 < NLO; n2++) {
				std::swap(X[n2 + NLO * n1], X[n2 + NLO * nr]);
			}
		}
		if (n1 != N1 - 1) {
			int k = N1 / R;
			nr++;
			while (Rm1 * k < nr) {
				nr -= Rm1 * k;
				k /= R;
			}
			nr += k - 1;
		}

	}

}

const std::vector<int>& get_digit_reversal(int R, int N) {
	static thread_local std::unordered_map<int, std::map<int, std::vector<int>>>values;
	auto iter = values[R].find(N);
	if( iter == values[R].end()) {
		std::vector<int> P(N);
		const int width = std::lround(std::log(N)/std::log(R));
		for( int i = 0; i < N; i++) {
			int j = i;
			int l = 0;
			for( int k = 0; k < width; k++) {
				l *= R;
				l += j % R;
				j /= R;
			}
			P[i] = l;
		}
		values[R][N] = std::move(P);
		iter = values[R].find(N);
	}
	return iter->second;
}
/*
void scramble_hi(double* X, int R, int NHI, int NLO) {
	const auto& P = get_digit_reversal(R, NHI);
	for (int i = 0; i < NHI; i++) {
		const int j = P[i];
		if (j > i) {
			double* Xi = X + NLO * i;
			double* Xj = X + NLO * j;
			for (int k = 0; k < NLO; k++) {
				std::swap(Xj[k], Xi[k]);
			}
		}
	}
}*/


void scramble(double* X, int R, int N) {
	int N1 = N;
	while (N1 * N1 > N) {
		N1 /= R;
	}
	const int N2 = N / (N1 * N1);
//	scramble_hi(X, R, N1, N1 * N2);
//	transpose(X, N1, N2);
//	scramble_hi(X, R, N1, N1 * N2);
}

void scramble(std::complex<double>* X, int R, int N) {
	int N1 = N;
	while ((size_t) N1 * (size_t) N1 > (size_t) N) {
		N1 /= R;
	}
	const int N2 = N / (N1 * N1);
	scramble_hi(X, R, N1, (N1 * N2) << 1);
	transpose(X, N1, N2);
	scramble_hi(X, R, N1, (N1 * N2) << 1);
}
