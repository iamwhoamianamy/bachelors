#include "matrix_mult.hpp"
#include "testing_helpers.hpp"
#include "integration.hpp"
#include "math.hpp"

namespace test
{
   void testMultiplication()
   {
      std::vector<std::vector<real>> matrix1 = {
         {1, 2, 3},
         {4, 5, 6},
         {7, 8, 9}
      };

      std::vector<std::vector<real>> matrix2 = {
         {5, 2, 8},
         {6, 4, 9},
         {3, 7, 1}
      };

      std::vector<real> vec1 = { 1, 1, 1 };
      std::vector<real> vec2 = { 2, 3, 4 };

      std::cout << math::mult(matrix1, vec1) << std::endl;
      std::cout << math::mult(vec1, vec2) << std::endl;
      std::cout << math::mult(matrix1, matrix2) << std::endl;
   }

   void compareWithMatrixMultiplication()
   {
      int n = 10;
      Vector3 childCenter(3, 1, 2);
      Vector3 parentCenter(5, 3, -1);
      Vector3 translation = childCenter - parentCenter;
      Torus torus = createTorus();
      auto bq = readBasisQuadratures();
      auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

      auto expansion = math::calcIntegralContribution(quadratures, n, childCenter);

      auto translated1 = Harmonics::translateWithComplex(expansion, translation);
      auto translated2 = translateWithMatrix(expansion, translation);
   }

   HarmonicSeries<Vector3> translateWithMatrix(
      const HarmonicSeries<Vector3>& expansion,
      const Vector3& translation)
   {
      auto regular = Harmonics::realToComplex(Harmonics::calcRegularSolidHarmonics(expansion.order(), translation));
      ComplexMatrix regularMatrix = regularToMatrix(regular);

      auto xComponent = Harmonics::realToComplex(Harmonics::separateX(expansion));
      auto yComponent = Harmonics::realToComplex(Harmonics::separateY(expansion));
      auto zComponent = Harmonics::realToComplex(Harmonics::separateZ(expansion));

      return Harmonics::createFormXYZ(
         Harmonics::complexToReal(ComplexHarmonicSeries(math::mult(regularMatrix, xComponent.data()))),
         Harmonics::complexToReal(ComplexHarmonicSeries(math::mult(regularMatrix, yComponent.data()))),
         Harmonics::complexToReal(ComplexHarmonicSeries(math::mult(regularMatrix, zComponent.data()))));
   }

   Matrix<Complex> regularToMatrix(
      const ComplexHarmonicSeries& regular)
   {
      ComplexMatrix res(
         regular.elemCount(),
         std::vector<Complex>(regular.elemCount()));

      for(int l = 0; l <= regular.order(); l++)
      {
         for(int m = -l; m <= l; m++)
         {
            for(int lambda = 0; lambda <= l; lambda++)
            {
               int dl = l - lambda;

               for(int mu = -lambda; mu <= lambda; mu++)
               {
                  int dm = m - mu;

                  if(-dl <= dm && dm <= dl)
                  {
                     res[l * l + l + m][dl * dl + dl + dm] = 
                        regular.getHarmonic(lambda * lambda + lambda + mu) * 
                        Harmonics::strangeFactor(m, mu);
                  }
               }
            }
         }
      }

      return res;
   }

}