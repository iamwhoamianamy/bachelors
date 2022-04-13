#include <vector>
#include "vector3.cuh"
#include "harmonics.hpp"

namespace test
{
   template <class T>
   using Matrix = std::vector<std::vector<T>>;

   typedef Matrix<std::complex<real>> ComplexMatrix;
   typedef Matrix<real> RealMatrix;

   void testMultiplication();
   void compareWithMatrixMultiplication();

   HarmonicSeries<Vector3> translateWithMatrix(
      const HarmonicSeries<Vector3>& expansion,
      const Vector3& translation);

   template <class T>
   std::vector<T> mult(
      const Matrix<T>& a,
      const std::vector<T>& b)
   {
      std::vector<T> res(b.size());

      for(size_t y = 0; y < a.size(); y++)
      {
         for(size_t x = 0; x < b.size(); x++)
         {
            res[y] += a[y][x] * b[x];
         }
      }

      return res;
   }

   ComplexMatrix regularToMatrix(
      const ComplexHarmonicSeries& regular);
}
