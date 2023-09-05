#include <fft86/transpose.hpp>
#include <fft86/vec.hpp>

#include <functional>

static void transpose(double* X, int xi, int yi, int D1, int N, int M) {
	if (xi <= yi && xi < N && yi < N) {
		if (M > 4 || xi + M > N || yi + M > N) {
			const int Mo2 = M >> 1;
			transpose(X, xi, yi, D1, N, Mo2);
			transpose(X, xi + Mo2, yi, D1, N, Mo2);
			transpose(X, xi, yi + Mo2, D1, N, Mo2);
			transpose(X, xi + Mo2, yi + Mo2, D1, N, Mo2);
		} else if (M == 4) {
			if (xi == yi) {
				auto* xptr = X + D1 * xi + yi;
				v4df X0, X1, X2, X3;
				X0.load(xptr + 0 * D1);
				X1.load(xptr + 1 * D1);
				X2.load(xptr + 2 * D1);
				X3.load(xptr + 3 * D1);
				transpose_4x4(X0, X1, X2, X3);
				X0.store(xptr + 0 * D1);
				X1.store(xptr + 1 * D1);
				X2.store(xptr + 2 * D1);
				X3.store(xptr + 3 * D1);
			} else {
				auto* xptr = X + D1 * xi + yi;
				auto* yptr = X + D1 * yi + xi;
				v4df X0, X1, X2, X3;
				v4df Y0, Y1, Y2, Y3;
				X0.load(xptr + 0 * D1);
				X1.load(xptr + 1 * D1);
				X2.load(xptr + 2 * D1);
				X3.load(xptr + 3 * D1);
				Y0.load(yptr + 0 * D1);
				Y1.load(yptr + 1 * D1);
				Y2.load(yptr + 2 * D1);
				Y3.load(yptr + 3 * D1);
				transpose_4x4(X0, X1, X2, X3);
				transpose_4x4(Y0, Y1, Y2, Y3);
				Y0.store(xptr + 0 * D1);
				Y1.store(xptr + 1 * D1);
				Y2.store(xptr + 2 * D1);
				Y3.store(xptr + 3 * D1);
				X0.store(yptr + 0 * D1);
				X1.store(yptr + 1 * D1);
				X2.store(yptr + 2 * D1);
				X3.store(yptr + 3 * D1);
			}
		} else if (M == 2) {
			if (xi == yi) {
				auto* xptr = X + D1 * xi + yi;
				v2df X0, X1;
				X0.load(xptr + 0 * D1);
				X1.load(xptr + 1 * D1);
				transpose_2x2(X0, X1);
				X0.store(xptr + 0 * D1);
				X1.store(xptr + 1 * D1);
			} else {
				auto* xptr = X + D1 * xi + yi;
				auto* yptr = X + D1 * yi + xi;
				v2df X0, X1;
				v2df Y0, Y1;
				X0.load(xptr + 0 * D1);
				X1.load(xptr + 1 * D1);
				Y0.load(yptr + 0 * D1);
				Y1.load(yptr + 1 * D1);
				transpose_2x2(X0, X1);
				transpose_2x2(Y0, Y1);
				X0.store(yptr + 0 * D1);
				X1.store(yptr + 1 * D1);
				Y0.store(xptr + 0 * D1);
				Y1.store(xptr + 1 * D1);
			}
		} else if (M == 1) {
			if (xi != yi) {
				auto* xptr = X + D1 * xi + yi;
				auto* yptr = X + D1 * yi + xi;
				std::swap(*xptr, *yptr);
			}
		}
	}
}

static void transpose(std::complex<double>* X, int xi, int yi, int D1, int N, int M) {
	if (xi <= yi && xi < N && yi < N) {
		if (M > 2 || xi + M > N || yi + M > N) {
			const int Mo2 = M >> 1;
			transpose(X, xi, yi, D1, N, Mo2);
			transpose(X, xi + Mo2, yi, D1, N, Mo2);
			transpose(X, xi, yi + Mo2, D1, N, Mo2);
			transpose(X, xi + Mo2, yi + Mo2, D1, N, Mo2);
		} else if (M == 2) {
			if (xi == yi) {
				auto* xptr = X + D1 * xi + yi;
				v4df X0, X1;
				X0.load(xptr + 0 * D1);
				X1.load(xptr + 1 * D1);
				transpose_2x2x2(X0, X1);
				X0.store(xptr + 0 * D1);
				X1.store(xptr + 1 * D1);
			} else {
				auto* xptr = X + D1 * xi + yi;
				auto* yptr = X + D1 * yi + xi;
				v4df X0, X1, Y0, Y1;
				X0.load(xptr + 0 * D1);
				X1.load(xptr + 1 * D1);
				Y0.load(yptr + 0 * D1);
				Y1.load(yptr + 1 * D1);
				transpose_2x2x2(X0, X1);
				transpose_2x2x2(Y0, Y1);
				X0.store(yptr + 0 * D1);
				X1.store(yptr + 1 * D1);
				Y0.store(xptr + 0 * D1);
				Y1.store(xptr + 1 * D1);
			}
		} else if (M == 1) {
			if (xi != yi) {
				v2df* xptr = (v2df*) (X + D1 * xi + yi);
				v2df* yptr = (v2df*) (X + D1 * yi + xi);
				std::swap(*xptr, *yptr);
			}
		}
	}
}

void transpose(double* X, int N, int NMID) {
	int M = 1 << (std::ilogb(N - 1) + 1);
	for (int imid = 0; imid < NMID; imid++) {
		transpose(X + imid * N, 0, 0, NMID * N, N, M);
	}
}

extern "C" void transpose_complex(std::complex<double>* X, int, int, int);
void transpose(std::complex<double>* X, int N, int NMID) {
	int M = 1 << (std::ilogb(N - 1) + 1);
	for (int imid = 0; imid < NMID; imid++) {
		transpose_complex(X + imid * N, NMID * N, N, M);
	//	transpose(X + imid * N, 0, 0, NMID * N, N, M);
	}
}
