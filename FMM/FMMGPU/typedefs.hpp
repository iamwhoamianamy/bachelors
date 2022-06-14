#pragma once
#include <vector>
//#include <thrust/complex.h>
#include <cuComplex.h>
#include "real.hpp"

template <class T>
using Matrix = std::vector<std::vector<T>>;

//typedef thrust::complex<real> Complex;

#ifdef REAL_IS_FLOAT
typedef cuComplex Complex;
#else
typedef cuDoubleComplex Complex;
#endif

typedef Matrix<Complex> ComplexMatrix;
typedef Matrix<real> RealMatrix;

typedef unsigned int uint;

#define __all__ __device__ __host__