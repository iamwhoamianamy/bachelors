#include <complex>
#include <limits>
#include "spherical_harmonics.hpp"

Factorials Harmonics::_factorials;

Harmonics::Harmonics(int n, const Vector3& point) :
   n(n)
{
   calcSphericalHarmonics(point);
}

const HarmonicSeries<real>& Harmonics::sphericalHarmonics() const
{
   return _sphericalHarmonics;
}

HarmonicSeries<real> Harmonics::calcRegularSolidHarmonics(size_t n, Vector3 point)
{
   return calcSolidHarmonics(n, point, true);
}

HarmonicSeries<real> Harmonics::calcIrregularSolidHarmonics(size_t n, Vector3 point)
{
   return calcSolidHarmonics(n, point, false);
}

void Harmonics::calcSphericalHarmonics(const Vector3& point)
{
   _sphericalHarmonics = HarmonicSeries<real>(n);
   fillWithLegendrePolynomials(point.z);
   fillWithLegendrePolynomialDerivatives(point.z);
   mirrorLegendrePolynomialDerivatives(point.z);
   addComplex(point.x, point.y);
}


void Harmonics::fillWithLegendrePolynomials(real z)
{
   _sphericalHarmonics.getHarmonic(0, 0) = 1;
   _sphericalHarmonics.getHarmonic(1, 0) = z;

   for(size_t l = 2; l < n; l++)
   {
      _sphericalHarmonics.getHarmonic(l, 0) = calcLegendrePolynomial(l, z);
   }
}

void Harmonics::fillWithLegendrePolynomialDerivatives(real z)
{
   _sphericalHarmonics.getHarmonic(0, 0) = 1;
   _sphericalHarmonics.getHarmonic(1, -1) = 1;
   _sphericalHarmonics.getHarmonic(1, 1) = 1;

   for(size_t l = 2; l < n; l++)
   {
      for(size_t m = 1; m <= l; m++)
      {
         _sphericalHarmonics.getHarmonic(l, m) = z * _sphericalHarmonics.getHarmonic(l - 1, m) +
            (l + m - 1) * _sphericalHarmonics.getHarmonic(l - 1, m - 1);
      }
   }
}

void Harmonics::mirrorLegendrePolynomialDerivatives(real z)
{
   for(int l = 2; l < n; l++)
   {
      for(int m = -l; m <= l; m++)
      {
         _sphericalHarmonics.getHarmonic(l, m) = _sphericalHarmonics.getHarmonic(l, -m);
      }
   }
}

real Harmonics::calcLegendrePolynomial(int l, real z)
{
   return ((2 * l - 1) * z * _sphericalHarmonics.getHarmonic(l - 1, 0) -
           (l - 1) * _sphericalHarmonics.getHarmonic(l - 2, 0)) / (l);
}

void Harmonics::addComplex(real x, real y)
{
   std::complex<real> ephi1m(sqrt(2.0));
   std::complex<real> mult(x, y);

   for(int m = 1; m < n; m++)
   {
      ephi1m *= mult;

      for(int l = m; l < n; l++)
      {
         _sphericalHarmonics.getHarmonic(l, -m) *= ephi1m.real();
         _sphericalHarmonics.getHarmonic(l, m) *= ephi1m.imag();
      }
   }
}

HarmonicSeries<real> Harmonics::calcSolidHarmonics(size_t n,
                                                   Vector3 point,
                                                   bool isRegular)
{
#if defined REAL_IS_DOUBLE
   real r = point.length() + std::numeric_limits<double>::epsilon();
#endif // !REAL_IS_DOUBLE
#if defined REAL_IS_FLOAT
   real r = point.length() + std::numeric_limits<float>::epsilon();
#endif // !REAL_IS_FLOAT

   point /= r;

   auto solidlHarmonics = Harmonics(n, point).sphericalHarmonics();

   real mult = isRegular ? r : 1 / r;
   real curr = isRegular ? 1 : mult;

   for(int l = 0; l < n; l++)
   {
      for(int m = 0; m <= l; m++)
      {
         solidlHarmonics.getHarmonic(l, -m) *= curr * (isRegular ? 1.0 / _factorials[l + m] : _factorials[l - m]);
      }

      for(int m = 1; m <= l; m++)
      {
         solidlHarmonics.getHarmonic(l, m) *= curr * (isRegular ? 1.0 / _factorials[l + m] : _factorials[l - m]);
      }

      curr *= -mult;
   }

   return solidlHarmonics;
}

