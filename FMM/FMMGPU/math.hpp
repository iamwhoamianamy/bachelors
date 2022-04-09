#pragma once
#include <vector>
#include "real.hpp"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

namespace math
{
   const real PI = 3.14159265359;
   const real mu0 = 1.2566370614e-6;
   const real SQRT_2 = sqrt(2.0);
   __device__ const real R_SQRT_2 = 0.70710678118;

   real calcFactorial(int n);
   real calcBinomial(int k, int n);

   template <class T>
   __all__ int sign(T val)
   {
      return (T(0) < val) - (val < T(0));
   }
}