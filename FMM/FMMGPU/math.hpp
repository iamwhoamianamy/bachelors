#pragma once
#include <vector>
#include "real.hpp"

namespace math
{
   const real PI = 3.14159265359;
   const real mu0 = 1.2566370614e-6;
   const real SQRT_2 = sqrt(2.0);
   const real R_SQRT_2 = 1.0 / sqrt(2.0);

   real calcFactorial(int n);
   real calcBinomial(int k, int n);

   template <class T>
   int sign(T val)
   {
      return (T(0) < val) - (val < T(0));
   }
}

