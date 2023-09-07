/*
 * vec.h
 *
 *  Created on: Aug 22, 2023
 *      Author: dmarce1
 */

#pragma once

#include <immintrin.h>

#include <complex>
#include <utility>

class v4df {
	__m256d v;
public:
	static constexpr int size() {
		return 4;
	}
	inline v4df() {
	}
	inline v4df(const v4df& other) {
		v = other.v;
	}
	inline v4df(double other) {
		v = _mm256_set1_pd(other);
	}
	inline v4df& operator=(double other) {
		v = _mm256_set1_pd(other);
		return *this;
	}
	inline v4df& operator=(const v4df& other) {
		v = other.v;
		return *this;
	}
	inline v4df& store(double* ptr) {
		_mm256_storeu_pd(ptr, v);
		return *this;
	}
	inline v4df& load(const double* ptr) {
		v = _mm256_loadu_pd(ptr);
		return *this;
	}
	inline v4df& store(std::complex<double>* ptr) {
		_mm256_storeu_pd((double*) ptr, v);
		return *this;
	}
	inline v4df& load(const std::complex<double>* ptr) {
		v = _mm256_loadu_pd((double*) ptr);
		return *this;
	}
	inline v4df& operator+=(const v4df& other) {
		v = _mm256_add_pd(v, other.v);
		return *this;
	}
	inline v4df& operator-=(const v4df& other) {
		v = _mm256_sub_pd(v, other.v);
		return *this;
	}
	inline v4df& operator*=(const v4df& other) {
		v = _mm256_mul_pd(v, other.v);
		return *this;
	}
	inline v4df& operator/=(const v4df& other) {
		v = _mm256_div_pd(v, other.v);
		return *this;
	}
	inline v4df operator+(const v4df& other) const {
		v4df rc = *this;
		rc += other;
		return rc;
	}
	inline v4df operator-(const v4df& other) const {
		v4df rc = *this;
		rc -= other;
		return rc;
	}
	inline v4df operator*(const v4df& other) const {
		v4df rc = *this;
		rc *= other;
		return rc;
	}
	inline v4df operator/(const v4df& other) const {
		v4df rc = *this;
		rc /= other;
		return rc;
	}
	inline v4df operator-() const {
		v4df rc = v4df(0.0) - *this;
		return rc;
	}
	inline double& operator[](int i) {
		return v[i];
	}
	inline double operator[](int i) const {
		return v[i];
	}
	friend inline void transpose_4x4(v4df& X0, v4df& X1, v4df& X2, v4df& X3) {
		v4df Y0, Y1, Y2, Y3;
		Y0.v = _mm256_unpacklo_pd(X0.v, X1.v);
		Y1.v = _mm256_unpackhi_pd(X0.v, X1.v);
		Y2.v = _mm256_unpacklo_pd(X2.v, X3.v);
		Y3.v = _mm256_unpackhi_pd(X2.v, X3.v);
		X0.v = _mm256_permute2f128_pd(Y0.v, Y2.v, 0x20);
		X1.v = _mm256_permute2f128_pd(Y1.v, Y3.v, 0x20);
		X2.v = _mm256_permute2f128_pd(Y0.v, Y2.v, 0x31);
		X3.v = _mm256_permute2f128_pd(Y1.v, Y3.v, 0x31);
	}

	friend inline void transpose_2x2x2(v4df& X0, v4df& X1) {
		v4df Y0;
		Y0.v = _mm256_permute2f128_pd(X0.v, X1.v, 0x20);
		X1.v = _mm256_permute2f128_pd(X0.v, X1.v, 0x31);
		X0 = Y0;
	}
	friend inline v4df complex_multiply(v4df a, v4df b) {
		v4df c;
		const auto res = _mm256_fmaddsub_pd(_mm256_shuffle_pd(a.v, a.v, 0), b.v,
				_mm256_mul_pd(_mm256_shuffle_pd(a.v, a.v, 15), _mm256_shuffle_pd(b.v, b.v, 5)));
		c.v = res;
		return c;
	}
	friend inline v4df rotate90(v4df a) {
		v4df b;
		v4df neg_odd;
		neg_odd.v = __m256d { 1.0, -1.0, 1.0, -1.0 };
		b.v = _mm256_permute4x64_pd(a.v, 177);
		b *= neg_odd;
		return b;
	}
};

class v2df {
	__m128d v;
public:
	static constexpr int size() {
		return 2;
	}
	inline v2df() {
	}
	inline v2df(const v2df& other) {
		v = other.v;
	}
	inline v2df(double other) {
		v = _mm_set1_pd(other);
	}
	inline v2df& operator=(double other) {
		v = _mm_set1_pd(other);
		return *this;
	}
	inline v2df& operator=(const v2df& other) {
		v = other.v;
		return *this;
	}
	inline v2df& store(double* ptr) {
		_mm_storeu_pd(ptr, v);
		return *this;
	}
	inline v2df& load(const double* ptr) {
		v = _mm_loadu_pd(ptr);
		return *this;
	}
	inline v2df& store(std::complex<double>* ptr) {
		_mm_storeu_pd((double*) ptr, v);
		return *this;
	}
	inline v2df& load(const std::complex<double>* ptr) {
		v = _mm_loadu_pd((double*) ptr);
		return *this;
	}
	inline v2df& operator+=(const v2df& other) {
		v = _mm_add_pd(v, other.v);
		return *this;
	}
	inline v2df& operator-=(const v2df& other) {
		v = _mm_sub_pd(v, other.v);
		return *this;
	}
	inline v2df& operator*=(const v2df& other) {
		v = _mm_mul_pd(v, other.v);
		return *this;
	}
	inline v2df& operator/=(const v2df& other) {
		v = _mm_div_pd(v, other.v);
		return *this;
	}
	inline v2df operator+(const v2df& other) const {
		v2df rc = *this;
		rc += other;
		return rc;
	}
	inline v2df operator-(const v2df& other) const {
		v2df rc = *this;
		rc -= other;
		return rc;
	}
	inline v2df operator*(const v2df& other) const {
		v2df rc = *this;
		rc *= other;
		return rc;
	}
	inline v2df operator/(const v2df& other) const {
		v2df rc = *this;
		rc /= other;
		return rc;
	}
	inline v2df operator-() const {
		v2df rc = v2df(0.0) - *this;
		return rc;
	}
	inline double& operator[](int i) {
		return v[i];
	}
	inline double operator[](int i) const {
		return v[i];
	}
	friend inline void transpose_2x2(v2df& X0, v2df& X1) {
		std::swap(X0.v[1], X1.v[0]);
	}
};

class v1df {
	double v;
public:
	static constexpr int size() {
		return 1;
	}
	inline v1df() {
	}
	inline v1df(const v1df& other) {
		v = other.v;
	}
	inline v1df(double other) {
		v = other;
	}
	inline v1df& operator=(double other) {
		v = other;
		return *this;
	}
	inline v1df& operator=(const v1df& other) {
		v = other.v;
		return *this;
	}
	inline v1df& store(double* ptr) {
		*ptr = v;
		return *this;
	}
	inline v1df& load(const double* ptr) {
		v = *ptr;
		return *this;
	}
	inline v1df& operator+=(const v1df& other) {
		v = (v + other.v);
		return *this;
	}
	inline v1df& operator-=(const v1df& other) {
		v = (v - other.v);
		return *this;
	}
	inline v1df& operator*=(const v1df& other) {
		v = (v * other.v);
		return *this;
	}
	inline v1df& operator/=(const v1df& other) {
		v = (v / other.v);
		return *this;
	}
	inline v1df operator+(const v1df& other) const {
		v1df rc = *this;
		rc += other;
		return rc;
	}
	inline v1df operator-(const v1df& other) const {
		v1df rc = *this;
		rc -= other;
		return rc;
	}
	inline v1df operator*(const v1df& other) const {
		v1df rc = *this;
		rc *= other;
		return rc;
	}
	inline v1df operator/(const v1df& other) const {
		v1df rc = *this;
		rc /= other;
		return rc;
	}
	inline v1df operator-() const {
		v1df rc(-v);
		return rc;
	}
	inline double& operator[](int i) {
		return v;
	}
	inline double operator[](int i) const {
		return v;
	}
	friend inline void transpose_2x2(v1df& X0, v1df& X1) {
	}
};
