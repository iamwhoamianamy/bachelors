#pragma once
#include <vector>
#include <complex>
#include "real.hpp"

template <class T>
using Matrix = std::vector<std::vector<T>>;

typedef Matrix<std::complex<real>> ComplexMatrix;
typedef Matrix<real> RealMatrix;