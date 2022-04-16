#pragma once
#include <vector>
#include <thrust/complex.h>
#include "real.hpp"

template <class T>
using Matrix = std::vector<std::vector<T>>;

typedef Matrix<thrust::complex<real>> ComplexMatrix;
typedef Matrix<real> RealMatrix;