#pragma once
#include <vector>
#include "complex.cuh"
#include "real.hpp"

template <class T>
using Matrix = std::vector<std::vector<T>>;

typedef Matrix<Complex> ComplexMatrix;
typedef Matrix<real> RealMatrix;