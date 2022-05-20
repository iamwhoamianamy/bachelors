#pragma once
#include <vector>
#include <thrust/complex.h>
#include "real.hpp"
#include "vector3.cuh"
#include "factorials.hpp"
#include "harmonic_series.hpp"

typedef HarmonicSeries<Complex> ComplexHarmonicSeries;
typedef HarmonicSeries<real> RealHarmonicSeries;

class Harmonics
{
private:
   size_t _order;
   RealHarmonicSeries _sphericalHarmonics;
   static Factorials _factorials;
public:   
   Harmonics(size_t order, const Vector3& point);

   const RealHarmonicSeries& sphericalHarmonics() const;

   static RealHarmonicSeries calcSolidHarmonics(
      size_t order, Vector3 point, bool isRegular);
   static RealHarmonicSeries calcRegularSolidHarmonics(
      size_t order, Vector3 point);
   static RealHarmonicSeries calcIrregularSolidHarmonics(
      size_t order, Vector3 point);

   static ComplexHarmonicSeries realToComplex(
      const RealHarmonicSeries& harmonics);

   static RealHarmonicSeries complexToReal(
      const ComplexHarmonicSeries& harmonics);

   static ComplexMatrix calcRealToComplexMatrix2D(size_t order);
   static ComplexMatrix calcComplexToRealMatrix2D(size_t order);

   static std::vector<Complex> calcRealToComplexTransitionMatrix1D(size_t order);
   static std::vector<Complex> calcComplexToRealTransitionMatrix1D(size_t order);

   static real getFactorial(size_t n);
   
   static RealHarmonicSeries separateCoord(
      const HarmonicSeries<Vector3>& harmonics,
      size_t i);

   static HarmonicSeries<Vector3> createFormXYZ(
      const RealHarmonicSeries& xs,
      const RealHarmonicSeries& ys,
      const RealHarmonicSeries& zs);

private:
   void calcSphericalHarmonics(const Vector3& point);
   void fillWithLegendrePolynomials(real z);
   void fillWithLegendrePolynomialDerivatives(real z);
   void mirrorLegendrePolynomialDerivatives(real z);
   real calcLegendrePolynomial(int i, real z);
   void addComplex(real x, real y);
};