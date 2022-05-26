#include "integration.hpp"
#include "harmonics.hpp"
#include "math.hpp"

namespace math
{
   real calcLegendrePolynomial(real x, int n)
   {
      switch(n)
      {
         case 0: return 1;
         case 1: return x;
         default: return ((2 * n - 1) * x * 
                          calcLegendrePolynomial(x, n - 1) + (1 - n) * 
                          calcLegendrePolynomial(x, n - 2)) / (n);
      }
   }

   Vector3 calcVectorFunctionIntegral(
      Vector3(*f)(const Vector3&, const Vector3&),
      const Vector3& point,
      const std::vector<Quadrature>& quadratures)
   {
      Vector3 res = 0;

      for(auto& quadrature : quadratures)
      {
         res += f(point, quadrature) * quadrature.weight;
      }

      return res;
   }

   Vector3 pointFromBasisQuadrature(
      const Tetrahedron& tetr,
      const Vector3& quadr)
   {
      return tetr.a() +
         (tetr.b() - tetr.a()) * quadr.x +
         (tetr.c() - tetr.a()) * quadr.y +
         (tetr.d() - tetr.a()) * quadr.z;
   }

   Vector3 calcAViaSimpleIntegration(
      real current,
      const Vector3& point,
      const std::vector<Tetrahedron>& mesh,
      const BasisQuadratures& basisQuadratures)
   {
      return calcVectorFunctionIntegral(
         simpleIntegrationFunctionForA,
         point, 
         tetrahedraToQuadratures(mesh, basisQuadratures)) / 
         (4 * PI) * current;
   }

   Vector3 simpleIntegrationFunctionForA(
      const Vector3& point,
      const Vector3& integr)
   {
      return integr.perp().normalized() / (point - integr).length();
   }

   Vector3 calcAViaMultipoleMethod(
      real current,
      const Vector3& point,
      const std::vector<Tetrahedron>& mesh,
      const BasisQuadratures& basisQuadratures,
      int harmonicOrder)
   {
      auto quadratures = tetrahedraToQuadratures(mesh, basisQuadratures);
      HarmonicSeries<Vector3> integrals = 
         calcIntegralContribution(quadratures, harmonicOrder);

      auto irregularHarmonics = 
         Harmonics::calcSolidHarmonics(harmonicOrder, point, false);

      Vector3 res;

      for(int l = 0; l <= harmonicOrder; l++)
      {
         for(int m = -l; m <= l; m++)
         {
            res += integrals.getHarmonic(l, m) * irregularHarmonics.getHarmonic(l, m);
         }
      }

      return res / (4.0 * PI) * current;
   }

   Vector3 calcBioSavartLaplace(
      real current,
      const Vector3& point,
      std::vector<Quadrature>& quadratures)
   {
      return calcVectorFunctionIntegral(
         bioSavartLaplaceFunction,
         point,
         quadratures) *
         MU0 / (4 * PI) * current;
   }

   Vector3 bioSavartLaplaceFunction(
      const Vector3& point,
      const Vector3& integr)
   {
      Vector3 diff = point - integr;
      return Vector3::cross(integr.perp().normalized(), (diff)) / 
         pow(diff.length(), 3);
   }

   HarmonicSeries<Vector3> calcIntegralContribution(
      std::vector<Quadrature>& quadratures,
      int harmonicOrder,
      const Vector3& center)
   {
      std::vector<Quadrature*> newQuadratures(quadratures.size());

      for(size_t i = 0; i < quadratures.size(); i++)
      {
         newQuadratures[i] = &quadratures[i];
      }

      return calcIntegralContribution(newQuadratures, harmonicOrder, center);
   }

   HarmonicSeries<Vector3> calcIntegralContribution(
      const std::vector<Quadrature*>& quadratures, 
      int harmonicOrder,
      const Vector3& center)
   {
      HarmonicSeries<Vector3> res(harmonicOrder);

      for(auto quadrature : quadratures)
      {
         auto regularHarmonics = Harmonics::calcSolidHarmonics(
            harmonicOrder,
            *quadrature - center,
            true);

         for(int l = 0; l <= harmonicOrder; l++)
         {
            for(int m = -l; m <= l; m++)
            {
               res.getHarmonic(l, m) +=
                  quadrature->perp().normalized() *
                  regularHarmonics.getHarmonic(l, m) * quadrature->weight;
            }
         }
      }

      return res;
   }

   std::vector<Quadrature> tetrahedraToQuadratures(
      const std::vector<Tetrahedron>& mesh,
      const BasisQuadratures& basisQuadratures)
   {
      auto res = std::vector<Quadrature>(mesh.size() * basisQuadratures.order());

      for(size_t t = 0; t < mesh.size(); t++)
      {
         real volume = mesh[t].volume();

         for(size_t i = 0; i < basisQuadratures.order(); i++)
         {
            Vector3 quadrature = math::pointFromBasisQuadrature(
               mesh[t],
               basisQuadratures.values(i));

            res[t * basisQuadratures.order() + i] = Quadrature(
               quadrature,
               volume,
               basisQuadratures.w(i));
         }
      }

      return res;
   }
}