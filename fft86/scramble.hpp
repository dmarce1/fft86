/*
 * scramble.hpp
 *
 *  Created on: Aug 22, 2023
 *      Author: dmarce1
 */

#ifndef SCRAMBLE_HPP_
#define SCRAMBLE_HPP_

#include <complex>
#include <vector>

#include <fft86/transpose.hpp>
#include <fft86/vec.hpp>


void scramble(std::complex<double>* X, int R, int N);
void scramble(double* X, int R, int N);
const std::vector<int>& get_digit_reversal(int R, int N);
//void scramble_hi(std::complex<double>* X, int R, int NHI, int NLO);
extern "C" void scramble_hi(std::complex<double>* X, int R, int NHI, int NLO);
void scramble_hi(double* X, int R, int N1, int NLO);


#endif /* SCRAMBLE_HPP_ */
