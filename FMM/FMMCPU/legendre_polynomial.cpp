#include "legendre_polynomial.hpp"

namespace legendre_polynomial
{
   real pn(real x, int n)
   {
      switch(n)
      {
         case 0: return 1;
         case 1: return x;
         default: return ((2 * n - 1) * x * pn(x, n - 1) + (1 - n)*pn(x, n - 2)) / (n);
      }
   }
}

