#include <iostream>
#include <complex>
#include "multipole_solver.hpp"
#include "math.hpp"
#include "spherical_harmonics.hpp"
#include "math.hpp"

MultipoleSolver::MultipoleSolver(const std::vector<Tetrahedron>& mesh, 
                                 const BasisQuadratures& basisQuadratures)
{
   _quadratures = math::tetrahedraToQuadratures(mesh, basisQuadratures);
   octreeRoot = new Octree(Box(Vector3(0, 0, 0), Vector3(3, 3, 3)), 10);
   octreeRoot->insert(_quadratures);
}

Vector3 MultipoleSolver::calcAFromRoot(real current, const Vector3& point)
{
   HarmonicSeries<Vector3> integralContribution = 
      math::calcIntegralContribution(_quadratures, _n);

   auto irregularHarmonics = Harmonics::calcSolidHarmonics(_n, point, false);
   Vector3 res;

   for(int l = 0; l < _n; l++)
   {
      for(int m = -l; m <= l; m++)
      {
         res += integralContribution.getHarmonic(l, m) *
            irregularHarmonics.getHarmonic(l, m);
      }
   }

   return res / (4.0 * math::PI) * current;
}

Vector3 MultipoleSolver::calcAWithMultipoleMethod(real current, const Vector3& point)
{
   octreeRoot->calcLocalMultipolesWithoutTranslation(_n);

   Vector3 res = octreeRoot->calcA(point);

   return res / (4.0 * math::PI) * current;
}

Vector3 MultipoleSolver::calcAwithFastMultipoleMethod(real current, const Vector3& point)
{
   for(auto child : octreeRoot->children())
   {
      child->calcLocalMultipolesWithoutTranslation(_n);
   }

   HarmonicSeries<Vector3> integralContribution(_n);

   //for(auto child : octreeRoot->children())
   //{
   //   Vector3 translation = child->box().center;

   //   for(int j = 0; j < _n; j++)
   //   {
   //      for(int k = -j; k <= j; k++)
   //      {
   //         Vector3 tempSum;
   //         auto R = Harmonics::calcSolidHarmonics(_n, translation, true);

   //         for(int l = 0; l <= j; l++)
   //         {
   //            for(int m = -l; m <= l; m++)
   //            {
   //               if(l - j <= k - m && k - m <= j - l)
   //               {
   //                  Vector3 o = child->_multipoleExpansion.getHarmonic(j - l, k - m);
   //                  //Vector3 o;
   //                  real i = std::pow(std::complex<real>(0, 1),
   //                                    abs(k) - abs(m) - abs(k - m)).real();
   //                  real a0 = math::calcAlm(l, m);
   //                  real a1 = math::calcAlm(j - l, k - m);
   //                  real a2 = math::calcAlm(j, k);

   //                  tempSum += o * i * a0 * a1 * R.getHarmonic(l, -m) / a2;
   //               }
   //            }
   //         }

   //         integralContribution.getHarmonic(j, k) += tempSum;
   //      }
   //   }
   //}

   for(auto child : octreeRoot->children())
   {
      integralContribution.add(child->multipoleExpansion());
   }

   auto irregularHarmonics = Harmonics::calcSolidHarmonics(_n, point, false);
   Vector3 res;

   for(int l = 0; l < _n; l++)
   {
      for(int m = -l; m <= l; m++)
      {
         res += integralContribution.getHarmonic(l, m) *
            irregularHarmonics.getHarmonic(l, m);
      }
   }

   return res / (4.0 * math::PI) * current;
}

MultipoleSolver::~MultipoleSolver()
{
   delete octreeRoot;
}
