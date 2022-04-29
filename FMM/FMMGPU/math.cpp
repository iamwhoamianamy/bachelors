#include "math.hpp"
#include "harmonics.hpp"
#include "iomanip"

namespace math
{
   real calcFactorial(int n)
   {
      return n <= 0 ? 1 : n * calcFactorial(n - 1);
   }

   real calcBinomial(int k, int n)
   {
       return Harmonics::getFactorial(n) / (Harmonics::getFactorial(k) * 
                                            Harmonics::getFactorial(n - k));
   }

   size_t nextDevisible(const size_t number, const size_t devidor)
   {
      if(devidor == 0)
         return number;
      else
         return (number + devidor - 1) / devidor * devidor;
   }
}

Complex& operator*(const Complex& lhs, const Complex& rhs)
{
   return cuCmulf(lhs, rhs);
}

__all__ Complex& operator*(const Complex& lhs, const real rhs)
{
   return make_cuComplex(cuCrealf(lhs) * rhs, cuCimagf(lhs) * rhs);
}

Complex& operator+(const Complex& lhs, const Complex& rhs)
{
   return cuCaddf(lhs, rhs);
}

__all__ Complex& operator-(const Complex& lhs, const Complex& rhs)
{
   return cuCsubf(lhs, rhs);
}

Complex& operator*=(Complex& lhs, const Complex& rhs)
{
   lhs = lhs * rhs;
   return lhs;
}

Complex& operator+=(Complex& lhs, const Complex& rhs)
{
   lhs = lhs + rhs;
   return lhs;
}

std::ostream& operator<<(std::ostream& os, const Complex& val)
{
   os << "(" << std::setw(10) << val.x << ", " << std::setw(10) << val.y << ")";
   return os;
}
