#include <complex>
#include <unordered_map>
#include <vector>
#include <memory>
#include <cassert>

const std::complex<double>* get_twiddles(int N) {
	static thread_local std::unordered_map<int, std::vector<std::complex<double>>>values;
	auto i = values.find(N);
	if (i == values.end()) {
		std::vector<std::complex<double>> W(N);
		const double c = -2.0 * M_PI / N;
		for (int n = 0; n < N; n++) {
			W[n] = std::polar(1.0, c * n);
		}
		values[N] = std::move(W);
		i = values.find(N);
	}
	return i->second.data();
}


const double* get_cos_twiddles(int N) {
	static thread_local std::unordered_map<int, std::vector<double>> values;
	auto i = values.find(N);
	if (i == values.end()) {
		std::vector<double> W(N);
		const double c = -2.0 * M_PI / N;
		for (int n = 0; n < N; n++) {
			W[n] = std::cos(c * n);
		}
		values[N] = std::move(W);
		i = values.find(N);
	}
	return i->second.data();
}

const double* get_sin_twiddles(int N) {
	static thread_local std::unordered_map<int, std::vector<double>> values;
	auto i = values.find(N);
	if (i == values.end()) {
		std::vector<double> W(N);
		const double c = -2.0 * M_PI / N;
		for (int n = 0; n < N; n++) {
			W[n] = std::sin(c * n);
		}
		values[N] = std::move(W);
		i = values.find(N);
	}
	return i->second.data();
}


const double* get_rotated_cos_twiddles(int R, int N) {
	static thread_local std::unordered_map<int,std::unordered_map<int, std::vector<double>>> values;
	auto i = values[R].find(N);
	if (i == values[R].end()) {
		std::vector<double> W(N);
		const double c = -2.0 * R * M_PI / N;
		for (int n = 0; n < N; n++) {
			W[n] = std::cos(c * n);
		}
		values[R][N] = std::move(W);
		i = values[R].find(N);
	}
	return i->second.data();
}

const double* get_rotated_sin_twiddles(int R, int N) {
	static thread_local std::unordered_map<int,std::unordered_map<int, std::vector<double>>> values;
	auto i = values[R].find(N);
	if (i == values[R].end()) {
		std::vector<double> W(N);
		const double c = -2.0 * R * M_PI / N;
		for (int n = 0; n < N; n++) {
			W[n] = std::sin(c * n);
		}
		values[R][N] = std::move(W);
		i = values[R].find(N);
	}
	return i->second.data();
}

std::shared_ptr<std::vector<double>> get_6step_cos_twiddles(int N1, int N2) {
	static thread_local std::unordered_map<int, std::unordered_map<int, std::shared_ptr<std::vector<double>>>> values;
	auto i = values[N1].find(N2);
	if (i == values[N1].end()) {
		std::vector<double> W(N1*N2);
		const double c = -2.0 * M_PI / (N1 * N2);
		for (int n1 = 0; n1 < N1; n1++) {
			for( int k2 = 0; k2 < N2; k2++) {
				W[n1+N1*k2] = std::cos(c * n1 * k2);
			}
		}
		values[N1][N2] = std::make_shared<std::vector<double>>(std::move(W));
		i = values[N1].find(N2);
		assert(i != values[N1].end());
	}
	return i->second;
}


std::shared_ptr<std::vector<double>> get_6step_sin_twiddles(int N1, int N2) {
	static thread_local std::unordered_map<int, std::unordered_map<int, std::shared_ptr<std::vector<double>>>> values;
	auto i = values[N1].find(N2);
	if (i == values[N1].end()) {
		std::vector<double> W(N1*N2);
		const double c = -2.0 * M_PI / (N1 * N2);
		for (int n1 = 0; n1 < N1; n1++) {
			for( int k2 = 0; k2 < N2; k2++) {
				W[n1+N1*k2] = std::sin(c * n1 * k2);
			}
		}
		values[N1][N2] = std::make_shared<std::vector<double>>(std::move(W));
		i = values[N1].find(N2);
		assert(i != values[N1].end());
	}
	return i->second;
}
