#include "math.hpp"
#include "harmonics.hpp"
#include "iomanip"

namespace math
{
   Vector3 cylindricToCartesian(const Vector3& point)
   {
      return  {
         point.y * std::cos(point.x),
         point.y* std::sin(point.x),
         point.z };
   }

   real calcFactorial(int n)
   {
      return n <= 0 ? 1 : n * calcFactorial(n - 1);
   }

   real calcBinomial(int k, int n)
   {
       return Harmonics::getFactorial(n) / (Harmonics::getFactorial(k) * 
                                            Harmonics::getFactorial(n - k));
   }

   real randBetween(real min, real max)
   {
      return (real)std::rand() / RAND_MAX * (max - min) + min;
   }

   size_t nextDevisible(const size_t number, const size_t devidor)
   {
      if(devidor == 0)
         return number;
      else
         return (number + devidor - 1) / devidor * devidor;
   }

   Box getBoundingBox(const std::vector<Vector3>& points)
   {
      real maxX = math::max(points, 0);
      real maxY = math::max(points, 1);
      real maxZ = math::max(points, 2);

      real minX = math::min(points, 0);
      real minY = math::min(points, 1);
      real minZ = math::min(points, 2);

      Vector3 halfDim(
         (maxX - minX) / 2 + 1e-5,
         (maxY - minY) / 2 + 1e-5,
         (maxZ - minZ) / 2 + 1e-5);

      Vector3 center(minX + halfDim.x, minY + halfDim.y, minZ + halfDim.z);

      return { center, halfDim * 1.1 };
   }
}

real math::max(const std::vector<Vector3>& vec, size_t axis)
{
   real res = -1e+16;

   for(size_t i = 0; i < vec.size(); i++)
   {
      if(vec[i][axis] > res)
         res = vec[i][axis];
   }

   return res;
}

real math::min(const std::vector<Vector3>& vec, size_t axis)
{
   real res = 1e+16;

   for(size_t i = 0; i < vec.size(); i++)
   {
      if(vec[i][axis] < res)
         res = vec[i][axis];
   }

   return res;
}

Complex operator*(const Complex& lhs, const Complex& rhs)
{
   return cuCmulf(lhs, rhs);
}

__all__ Complex operator*(const Complex& lhs, const real rhs)
{
   return make_cuComplex(cuCrealf(lhs) * rhs, cuCimagf(lhs) * rhs);
}

Complex operator+(const Complex& lhs, const Complex& rhs)
{
   return cuCaddf(lhs, rhs);
}

__all__ Complex operator-(const Complex& lhs, const Complex& rhs)
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
