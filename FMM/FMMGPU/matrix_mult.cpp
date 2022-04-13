#include "matrix_mult.hpp"

void compareWithMatrixMultiplication()
{
   int n = 10;
   Vector3 childCenter(3, 1, 2);
   Vector3 parentCenter(6, 4, -1);

   Vector3 translation = childCenter - parentCenter;

   auto regular = Harmonics::calcRegularSolidHarmonics(n, translation);


}