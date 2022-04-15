#include <vector>
#include "vector3.cuh"
#include "harmonics.hpp"
#include "typedefs.hpp"

namespace test
{
   void testMultiplication();
   void testBLASMultiplication();
   void compareWithMatrixMultiplication();

   HarmonicSeries<Vector3> translateWithMatrix(
      const HarmonicSeries<Vector3>& expansion,
      const Vector3& translation);

   ComplexMatrix regularToMatrix(
      const ComplexHarmonicSeries& regular);
}
