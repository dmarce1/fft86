#include <utility>
#include <vector>

template<class T>
void transpose_ruritarian(T* X, int N1, int N2) {
	const int N = N1 * N2;
	if (N1 > N2) {
		std::swap(N1, N2);
	}

	for (int n1 = 0; n1 < N1; n1++) {
		std::vector<bool> visited(false, N2);
		for (int n2 = 0; n2 < N2; n2++) {
			if (!visited[n2]) {
				int p = n2;
				do {
					const int q = (n2 * N1) % N2;
					visited[q] = true;
					if (q != n2) {
						std::swap(X[n1 * N2 + p], X[n1 * N2 + q]);
					}
					q = p;
				} while (p != n2);
			}
		}
	}
	for (int O = 0; O < N; O += N1 * N1) {
		for (int n1 = 0; n1 < N1; n1++) {
			for (int n2 = 0; n2 < N1; n2++) {
				std::swap(X[N1 * n2 + n1], X[N1 * n1 + n2]);
			}
		}

		for (int n1 = 0; n1 < N1; n1++) {
			for (int n2 = 0; n2 < N1; n2++) {
				std::swap(X[N1 * n2 + n1], X[N1 * n1 + n2]);
			}
		}
	}
}
