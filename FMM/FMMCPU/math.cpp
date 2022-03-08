#include "math.hpp"

namespace math
{
   real calcLegendrePolynomial(real x, int n)
   {
      switch(n)
      {
         case 0: return 1;
         case 1: return x;
         default: return ((2 * n - 1) * x * calcLegendrePolynomial(x, n - 1) + (1 - n)*calcLegendrePolynomial(x, n - 2)) / (n);
      }
   }

   real calcFactorial(int n)
   {
      return n <= 0 ? 1 : n * calcFactorial(n - 1);
   }
}

