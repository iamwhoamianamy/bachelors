#pragma once
#include <vector>
#include "real.hpp"
#include "vector3.hpp"
#include "factorials.hpp"
#include "harmonic_series.hpp"

class Harmonics
{
private:
   int _n;
   HarmonicSeries<real> _sphericalHarmonics;
   static Factorials _factorials;
public:   
   Harmonics(int n, const Vector3& point);

   const HarmonicSeries<real>& sphericalHarmonics() const;

   static HarmonicSeries<real> calcRegularSolidHarmonics(size_t n,
                                                         Vector3 point);
   static HarmonicSeries<real> calcIrregularSolidHarmonics(size_t n,
                                                           Vector3 point);
   static HarmonicSeries<real> calcSolidHarmonics(size_t n,
                                                  Vector3 point,
                                                  bool isRegular);
   template <class T>
   static T getHarmonic(int l, int m, const std::vector<std::vector<T>>& harmonics);

   static real getFactorial(size_t n)
   {
      return _factorials[n];
   }

private:
   void calcSphericalHarmonics(const Vector3& point);
   void fillWithLegendrePolynomials(real z);
   void fillWithLegendrePolynomialDerivatives(real z);
   void mirrorLegendrePolynomialDerivatives(real z);
   real calcLegendrePolynomial(int i, real z);
   void addComplex(real x, real y);
};

template <class T>
T Harmonics::getHarmonic(int l, int m, const std::vector<std::vector<T>>& harmonics)
{
   return harmonics[l][harmonics[l].size() / 2 + m];
}
