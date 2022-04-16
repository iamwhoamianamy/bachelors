#pragma once
#include <vector>
#include <thrust/complex.h>
#include "real.hpp"
#include "vector3.cuh"
#include "factorials.hpp"
#include "harmonic_series.hpp"

typedef HarmonicSeries<thrust::complex<real>> ComplexHarmonicSeries;
typedef HarmonicSeries<real> RealHarmonicSeries;

class Harmonics
{
private:
   size_t _order;
   RealHarmonicSeries _sphericalHarmonics;
   static Factorials _factorials;
public:   
   Harmonics(int order, const Vector3& point);

   const RealHarmonicSeries& sphericalHarmonics() const;

   static RealHarmonicSeries calcSolidHarmonics(
      size_t order,Vector3 point, bool isRegular);
   static RealHarmonicSeries calcRegularSolidHarmonics(
      size_t order, Vector3 point);
   static RealHarmonicSeries calcIrregularSolidHarmonics(
      size_t order, Vector3 point);

   static ComplexHarmonicSeries realToComplex(
      const RealHarmonicSeries& harmonics);
   static RealHarmonicSeries complexToReal(
      const ComplexHarmonicSeries& harmonics);

   static ComplexHarmonicSeries translate(
      const ComplexHarmonicSeries& a,
      const ComplexHarmonicSeries& b);

   static RealHarmonicSeries translate(
      const RealHarmonicSeries& a,
      const RealHarmonicSeries& b);

   static HarmonicSeries<Vector3> translateWithComplex(
      const HarmonicSeries<Vector3>& expansion,
      const Vector3& translation);

   static HarmonicSeries<Vector3> translateWithReal(
      const HarmonicSeries<Vector3>& expansion,
      const Vector3& translation);

   static real getFactorial(size_t n);
   
   static RealHarmonicSeries separateCoord(
      const HarmonicSeries<Vector3>& harmonics,
      size_t i);

   static RealHarmonicSeries separateX(
      const HarmonicSeries<Vector3>& harmonics);

   static RealHarmonicSeries separateY(
      const HarmonicSeries<Vector3>& harmonics);

   static RealHarmonicSeries separateZ(
      const HarmonicSeries<Vector3>& harmonics);

   static HarmonicSeries<Vector3> createFormXYZ(
      const RealHarmonicSeries& xs,
      const RealHarmonicSeries& ys,
      const RealHarmonicSeries& zs);

   static real strangeFactor(int m, int mu);
private:
   void calcSphericalHarmonics(const Vector3& point);
   void fillWithLegendrePolynomials(real z);
   void fillWithLegendrePolynomialDerivatives(real z);
   void mirrorLegendrePolynomialDerivatives(real z);
   real calcLegendrePolynomial(int i, real z);
   void addComplex(real x, real y);
};