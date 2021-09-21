#include "point.h"
#include "exeptions.h"

namespace triangle_quadratures
{
   double& Point::operator[](int i)
   {
      switch(i)
      {
         case 0: return x;
         case 1: return y;
         case 2: return z;
         default: throw RangeExeption();
      }
   }
}