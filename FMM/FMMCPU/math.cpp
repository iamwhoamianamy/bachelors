#include "math.hpp"
#include "spherical_harmonics.hpp"

namespace math
{
   real calcLegendrePolynomial(real x, int n)
   {
      switch(n)
      {
         case 0: return 1;
         case 1: return x;
         default: return ((2 * n - 1) * x * calcLegendrePolynomial(x, n - 1) + (1 - n)*calcLegendrePolynomial(x, n - 2)) / (n);
      }
   }

   real calcFactorial(int n)
   {
      return n <= 0 ? 1 : n * calcFactorial(n - 1);
   }

   Vector3 calcVectorFunctionIntegral(Vector3(*f)(const Vector3&, const Vector3&),
                                      const Vector3& point,
                                      const std::vector<Tetrahedron>& mesh,
                                      const BasisQuadratures& basisQuadratures)
   {
      Vector3 res = 0;

      for(auto& tetrahedron : mesh)
      {
         Vector3 tetrahedronRes = 0;

         for(size_t i = 0; i < basisQuadratures.order(); i++)
         {
            Vector3 quadrature = pointFromBasisQuadrature(tetrahedron, basisQuadratures.values(i));
            tetrahedronRes += f(point, quadrature) * basisQuadratures.w(i);
         }

         res += tetrahedronRes * tetrahedron.volume();
      }

      return res;
   }

   Vector3 pointFromBasisQuadrature(const Tetrahedron& tetr,
                                    const Vector3& quadr)
   {
      return tetr.a() +
         (tetr.b() - tetr.a()) * quadr.x +
         (tetr.c() - tetr.a()) * quadr.y +
         (tetr.d() - tetr.a()) * quadr.z;
   }

   Vector3 calcAViaSimpleIntegration(real current, const Vector3& point, 
                                     const std::vector<Tetrahedron>& mesh, 
                                     const BasisQuadratures& basisQuadratures)
   {
      return calcVectorFunctionIntegral(simpleIntegrationFunctionForA,
                                        point, mesh, basisQuadratures) / (4 * PI) * current;
   }

   Vector3 simpleIntegrationFunctionForA(const Vector3& point, const Vector3& integr)
   {
      return integr.perp().normalized() / (point - integr).length();
   }

   Vector3 calcAViaMultipoleMethod(real current, const Vector3& point,
                                   const std::vector<Tetrahedron>& mesh,
                                   const BasisQuadratures& basisQuadratures, int n)
   {
      auto quadratures = tetrahedraToQuadratures(mesh, basisQuadratures);
      HarmonicSeries<Vector3> integrals = calcIntegralContribution(quadratures, n);

      auto irregularHarmonics = Harmonics::calcSolidHarmonics(n, point, false);
      Vector3 res;

      for(int l = 0; l < n; l++)
      {
         for(int m = -l; m <= l; m++)
         {
            res += integrals[l][integrals[l].size() / 2 + m] *
               irregularHarmonics.getHarmonic(l, m);
         }
      }

      return res / (4.0 * PI) * current;
   }

   Vector3 calcBioSavartLaplace(real current, const Vector3& point,
                                const std::vector<Tetrahedron>& mesh,
                                const BasisQuadratures& basisQuadratures)
   {
      return calcVectorFunctionIntegral(bioSavartLaplaceFunction,
                                        point, mesh, basisQuadratures) *
         mu0 / (4 * PI) * current;
   }

   Vector3 bioSavartLaplaceFunction(const Vector3& point, const Vector3& integr)
   {
      Vector3 diff = point - integr;
      return Vector3::cross(integr.perp().normalized(), (diff)) / pow(diff.length(), 3);
   }

   HarmonicSeries<Vector3> calcIntegralContribution(std::vector<Quadrature>& quadratures,
                                                    int n, const Vector3& center)
   {
      std::vector<Quadrature*> newQuadratures(quadratures.size());

      for(size_t i = 0; i < quadratures.size(); i++)
      {
         newQuadratures[i] = &quadratures[i];
      }

      return calcIntegralContribution(newQuadratures, n, center);
   }

   HarmonicSeries<Vector3> calcIntegralContribution(const std::vector<Quadrature*>& quadratures,
                                                    int n, const Vector3& center)
   {
      HarmonicSeries<Vector3> res(n);

      for(auto quadrature : quadratures)
      {
         auto regularHarmonics = Harmonics::calcSolidHarmonics(n,
                                                               quadrature->coordinates - center,
                                                               true);
         for(int l = 0; l < n; l++)
         {
            for(int m = -l; m <= l; m++)
            {
               res.getHarmonic(l, m) +=
                  quadrature->coordinates.perp().normalized() *
                  regularHarmonics.getHarmonic(l, m) * quadrature->weight;
            }
         }
      }

      return res;
   }

   real calcAlm(int l, int m)
   {
      return pow(-1, l) / sqrt(Harmonics::getFactorial(l - m) *
                               Harmonics::getFactorial(l + m));
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
            Vector3 quadrature = math::pointFromBasisQuadrature(mesh[t],
                                                                basisQuadratures.values(i));
            res[t * basisQuadratures.order() + i] = Quadrature(quadrature, volume,
                                                               basisQuadratures.w(i));
         }
      }

      return res;
   }
}

