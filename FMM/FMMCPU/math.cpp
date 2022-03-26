#include "math.hpp"

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
            Vector3 quadrature = pointFromBasisQuadrature(tetrahedron, basisQuadratures.value()[0]);
            tetrahedronRes += f(point, quadrature) * basisQuadratures.w()[i];
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

   Vector3 calcAViaSimpleIntegration(real current, const Vector3& point, const std::vector<Tetrahedron>& mesh, const BasisQuadratures& basisQuadratures)
   {
      return calcVectorFunctionIntegral(simpleIntegrationFunctionForA,
                                        point, mesh, basisQuadratures) / (4 * PI) * current;
   }

   Vector3 simpleIntegrationFunctionForA(const Vector3& point, const Vector3& integr)
   {
      return 1 / (point - integr).length();
   }

   Vector3 calcAViaMultipoleMethod(real current, const Vector3& point, const std::vector<Tetrahedron>& mesh)
   {
      Vector3 res;



      return res;
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
}

