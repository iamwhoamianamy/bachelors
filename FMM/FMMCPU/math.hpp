#pragma once
#include <vector>
#include "real.hpp"
#include "vector3.hpp"
#include "tetrahedron.hpp"
#include "basis_quadratures.hpp"
#include "quadrature.hpp"
#include "harmonic_series.hpp"

namespace math
{
   const real PI = 3.14159265359;
   const real mu0 = 1.2566370614e-6;
   real calcLegendrePolynomial(real x, int n);
   real calcFactorial(int n);
   Vector3 calcVectorFunctionIntegral(Vector3 (*f)(const Vector3&, const Vector3&), 
                                      const Vector3& point,
                                      const std::vector<Quadrature>& quadratures);

   Vector3 pointFromBasisQuadrature(const Tetrahedron& tetr,
                                    const Vector3& quadr);

   Vector3 calcAViaSimpleIntegration(real current, const Vector3& point,
                                     const std::vector<Tetrahedron>& mesh,
                                     const BasisQuadratures& basisQuadratures);
   Vector3 simpleIntegrationFunctionForA(const Vector3& point,
                                         const Vector3& integr);

   Vector3 calcAViaMultipoleMethod(real current, const Vector3& point,
                                   const std::vector<Tetrahedron>& mesh,
                                   const BasisQuadratures& basisQuadratures, int n);

   Vector3 calcBioSavartLaplace(real current, const Vector3& point, std::vector<Quadrature>& quadratures);
   Vector3 bioSavartLaplaceFunction(const Vector3& point,
                                    const Vector3& integr);
   
   HarmonicSeries<Vector3> calcIntegralContribution(
      std::vector<Quadrature>& quadratures,
      int n, const Vector3& center = Vector3());

   HarmonicSeries<Vector3> calcIntegralContribution(
      const std::vector<Quadrature*>& quadratures,
      int n, const Vector3& center = Vector3());

   real calcAlm(int l, int m);

   std::vector<Quadrature> tetrahedraToQuadratures(
      const std::vector<Tetrahedron>& mesh,
      const BasisQuadratures& basisQuadratures);
}

