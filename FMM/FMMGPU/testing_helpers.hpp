#pragma once
#include <chrono>
#include <iostream>
#include "torus.hpp"
#include "basis_quadratures.hpp"
#include "exeptions.hpp"
#include "typedefs.hpp"

namespace test
{
   Torus createTorus();
   double getTime(void (*f)());
   double getTime(const std::chrono::steady_clock::time_point& start,
                  const std::chrono::steady_clock::time_point& stop);
   BasisQuadratures readBasisQuadratures();
   std::vector<Vector3> createPoints(const Vector3& begin, const Vector3& end, int n);
   
   template <class T>
   std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
   {
      for(size_t i = 0; i < vec.size(); i++)
      {
         os << vec[i] << " ";
      }

      return os;
   }

   template <class T>
   std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix)
   {
      for(size_t i = 0; i < matrix.size(); i++)
      {
         os << matrix[i] << std::endl;
      }

      return os;
   }
}
