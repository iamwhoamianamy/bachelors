#pragma once
#include <vector>
#include <complex>
#include "real.hpp"
#include "vector3.hpp"
#include "factorials.hpp"
#include "harmonic_series.hpp"

class Harmonics
{
private:
   int n;
   HarmonicSeries<real> _sphericalHarmonics;
   static Factorials _factorials;
public:   
   Harmonics(int n, const Vector3& point);

   const HarmonicSeries<real>& sphericalHarmonics() const;

   static HarmonicSeries<real> calcSolidHarmonics(
      size_t n,Vector3 point, bool isRegular);
   static HarmonicSeries<real> calcRegularSolidHarmonics(
      size_t n, Vector3 point);
   static HarmonicSeries<real> calcIrregularSolidHarmonics(
      size_t n, Vector3 point);

   static HarmonicSeries<std::complex<real>> realToComplex(
      const HarmonicSeries<real>& harmonics);
   static HarmonicSeries<real> complexToReal(
      const HarmonicSeries<std::complex<real>>& harmonics);

   static HarmonicSeries<std::complex<real>> translate(
      const HarmonicSeries<std::complex<real>>& a,
      const HarmonicSeries<std::complex<real>>& b);

   static HarmonicSeries<real> translate(
      const HarmonicSeries<real>& a,
      const HarmonicSeries<real>& b);

   static HarmonicSeries<Vector3> translateWithComplex(
      const HarmonicSeries<Vector3>& expansion,
      const Vector3& translation);

   static HarmonicSeries<Vector3> translateWithReal(
      const HarmonicSeries<Vector3>& expansion,
      const Vector3& translation);

   static real getFactorial(size_t n);
   
   static HarmonicSeries<real> separateX(const HarmonicSeries<Vector3>& harmonics);
   static HarmonicSeries<real> separateY(const HarmonicSeries<Vector3>& harmonics);
   static HarmonicSeries<real> separateZ(const HarmonicSeries<Vector3>& harmonics);
   static HarmonicSeries<Vector3> createFormXYZ(const HarmonicSeries<real>& xs,
                                                const HarmonicSeries<real>& ys,
                                                const HarmonicSeries<real>& zs);

   static real strangeFactor(int m, int mu);
private:
   void calcSphericalHarmonics(const Vector3& point);
   void fillWithLegendrePolynomials(real z);
   void fillWithLegendrePolynomialDerivatives(real z);
   void mirrorLegendrePolynomialDerivatives(real z);
   real calcLegendrePolynomial(int i, real z);
   void addComplex(real x, real y);
};
