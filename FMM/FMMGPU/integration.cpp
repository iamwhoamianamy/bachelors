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
      const Tetrahedron& tetrahedron,
      const Vector3& quadr)
   {
      return tetrahedron.a() +
         (tetrahedron.b() - tetrahedron.a()) * quadr.x +
         (tetrahedron.c() - tetrahedron.a()) * quadr.y +
         (tetrahedron.d() - tetrahedron.a()) * quadr.z;
   }

   Vector3 pointFromBasisQuadrature(
      const Triangle& triangle,
      const Vector3& quadr)
   {
      return triangle.a() +
         (triangle.b() - triangle.a()) * quadr.x +
         (triangle.c() - triangle.a()) * quadr.y;
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
         auto translation = *quadrature - center;
         auto regularHarmonics = Harmonics::calcSolidHarmonics(
            harmonicOrder,
            translation,
            true);

         for(size_t i = 0; i < regularHarmonics.elemCount(); i++)
         {
            res.getHarmonic(i) +=
               quadrature->perp().normalized() *
               regularHarmonics.getHarmonic(i) * quadrature->weight;
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

#pragma region BEM problem

   Vector3 BEMFunction(
      const Vector3& point,
      const BEMQuadrature& quadrature)
   {
      Vector3 diff = quadrature - point;
      real length3 = std::pow<real>(diff.length(), 3);
      Vector3 first = diff * Vector3::dot(quadrature.B, quadrature.normal);
      Vector3 second = Vector3::cross(diff,Vector3::cross(quadrature.B, quadrature.normal));

      return (first + second) / length3;
   }

   Vector3 calcBEMIntegral(
      const Vector3& point,
      const std::vector<BEMQuadrature>& quadratures)
   {
      Vector3 res = 0;

      for(auto& quadrature : quadratures)
      {
         res += BEMFunction(point, quadrature) * quadrature.weight;
      }

      return res / (4 * PI);
   }

   std::vector<BEMQuadrature> calcBEMquadraturesFromTriangles(
      const std::vector<Triangle>& triangles,
      const BasisQuadratures& basisQuadratures,
      const std::vector<ReferenceCylinderData>& referenceCylinderData,
      int normalDir)
   {
      std::vector<BEMQuadrature> result;
      result.reserve(basisQuadratures.order() * triangles.size());

      for (auto triangle : triangles)
      {
         Triangle triangleCartesian = Triangle(
            math::cylindricToCartesian(triangle.a()),
            math::cylindricToCartesian(triangle.b()),
            math::cylindricToCartesian(triangle.c()));

         //real area = triangle.area();
         real square = triangleCartesian.area();
         Vector3 B;
         Vector3 closestB;
         real minDist = 1000000;

         Box triangleBoundingBox = triangleCartesian.boundingBox();

         for (const auto & cylinderData : referenceCylinderData)
         {
            if(triangleBoundingBox.contains(cylinderData.point))
            {
               B = cylinderData.B;
               break;
            }
            else
            {
               real dist;

               if((dist = Vector3::distanceSquared(cylinderData.point, triangleCartesian.center())) < 
                  minDist * minDist)
               {
                  closestB = cylinderData.B;
                  minDist = dist;
               }
            }
         }

         if(B.x == 0 && B.y == 0 && B.z == 0)
            B = closestB;

         for(size_t q = 0; q < basisQuadratures.order(); q++)
         {
            auto pointFromBasisQuadratureCylindric = math::pointFromBasisQuadrature(
               triangle, basisQuadratures.values(q));

            Vector3 pointFromBasisQuadrature = 
               math::cylindricToCartesian(pointFromBasisQuadratureCylindric);

            Vector3 normal;

            if(normalDir == 0)
            {
               normal = pointFromBasisQuadrature;
               normal.z = 0;
               normal.normalize();
            }
            else
            {
               normal = Vector3::zAxis() * normalDir;
            }

            result.emplace_back(
               pointFromBasisQuadrature,
               square,
               basisQuadratures.w(q),
               B,
               normal);
         }
      }

      return result;
   }

   HarmonicSeries<Vector3> calcBEMIntegralContribution(
      const std::vector<Quadrature*>& quadratures,
      int harmonicOrder,
      const Vector3& center)
   {
      HarmonicSeries<Vector3> res(harmonicOrder);

      for(auto quadrature : quadratures)
      {
         auto translation = *quadrature - center;

         auto hx1 = Harmonics::calcRegularSolidHarmonics(
            harmonicOrder,
            translation + Vector3::xAxis() * math::eps);

         auto hx2 = Harmonics::calcRegularSolidHarmonics(
            harmonicOrder,
            translation - Vector3::xAxis() * math::eps);

         auto hy1 = Harmonics::calcRegularSolidHarmonics(
            harmonicOrder,
            translation + Vector3::yAxis() * math::eps);

         auto hy2 = Harmonics::calcRegularSolidHarmonics(
            harmonicOrder,
            translation - Vector3::yAxis() * math::eps);

         auto hz1 = Harmonics::calcRegularSolidHarmonics(
            harmonicOrder,
            translation + Vector3::zAxis() * math::eps);

         auto hz2 = Harmonics::calcRegularSolidHarmonics(
            harmonicOrder,
            translation - Vector3::zAxis() * math::eps);

         hx1.subtract(hx2);
         hy1.subtract(hy2);
         hz1.subtract(hz2);

         auto bemQuadrature = reinterpret_cast<BEMQuadrature*>(quadrature);

         for(size_t i = 0; i < hx1.elemCount(); i++)
         {
            Vector3 gradR = Vector3(
               hx1.getHarmonic(i),
               hy1.getHarmonic(i),
               hz1.getHarmonic(i)) / (2 * math::eps);

            Vector3 first = gradR * Vector3::dot(bemQuadrature->B, bemQuadrature->normal);
            Vector3 second = Vector3::cross(
               gradR,
               Vector3::cross(bemQuadrature->B, bemQuadrature->normal));

            res.getHarmonic(i) -= (first + second) * bemQuadrature->weight;
         }
      }

      return res;
   }

   real calcSurfaceArea(const std::vector<BEMQuadrature>& quadratures)
   {
      real result = 0;

      for(auto& quadrature : quadratures)
      {
         result += quadrature.weight;
      }

      return result;
   }

#pragma endregion
}