#pragma once
#include <vector>
#include <thrust/complex.h>
#include "real.hpp"

template <class T>
using Matrix = std::vector<std::vector<T>>;

typedef thrust::complex<real> Complex;

typedef Matrix<Complex> ComplexMatrix;
typedef Matrix<real> RealMatrix;

typedef unsigned int uint;