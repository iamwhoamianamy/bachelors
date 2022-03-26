#include <complex>
#include <limits>
#include "spherical_harmonics.hpp"

Factorials Harmonics::_factorials;

Harmonics::Harmonics(int n, const Vector3& point) :
   _n(n)
{
   calcSphericalHarmonics(point);
}

const std::vector<std::vector<real>>& Harmonics::sphericalHarmonics() const
{
   return _sphericalHarmonics;
}

void Harmonics::calcSphericalHarmonics(const Vector3& point)
{
   initHarmonicArrays();
   fillWithLegendrePolynomials(point.z);
   fillWithLegendrePolynomialDerivatives(point.z);
   mirrorLegendrePolynomialDerivatives(point.z);
   addComplex(point.x, point.y);
}

void Harmonics::initHarmonicArrays()
{
   _sphericalHarmonics = std::vector<std::vector<real>>(_n);

   for(size_t i = 0; i < _n; i++)
   {
      _sphericalHarmonics[i] = std::vector<real>(2 * i + 1);
   }
}

void Harmonics::fillWithLegendrePolynomials(real z)
{
   _sphericalHarmonics[0][0] = 1;
   _sphericalHarmonics[1][1] = z;

   for(size_t i = 2; i < _n; i++)
   {
      _sphericalHarmonics[i][i] = calcLegendrePolynomial(i, z);
   }
}

void Harmonics::fillWithLegendrePolynomialDerivatives(real z)
{
   _sphericalHarmonics[0][0] = 1;
   _sphericalHarmonics[1][0] = 1;
   _sphericalHarmonics[1][2] = 1;

   for(size_t l = 2; l < _n; l++)
   {
      for(size_t m = 1; m < l; m++)
      {
         size_t center = _sphericalHarmonics[l].size() / 2;
         size_t prevCenter = _sphericalHarmonics[l - 1].size() / 2;
         _sphericalHarmonics[l][center + m] = z * _sphericalHarmonics[l - 1][prevCenter + m] +
            (l + m - 1) * _sphericalHarmonics[l - 1][prevCenter + m - 1];
      }

      _sphericalHarmonics[l][_sphericalHarmonics[l].size() - 1] = (2 * l - 1) *
         _sphericalHarmonics[l - 1][_sphericalHarmonics[l - 1].size() - 1];
   }
}

void Harmonics::mirrorLegendrePolynomialDerivatives(real z)
{
   for(size_t l = 2; l < _n; l++)
   {
      for(size_t m = 0; m < _sphericalHarmonics[l].size() / 2; m++)
      {
         _sphericalHarmonics[l][m] = _sphericalHarmonics[l][_sphericalHarmonics[l].size() - 1 - m];
      }
   }
}

real Harmonics::calcLegendrePolynomial(int i, real z)
{
   return ((2 * i - 1) * z * _sphericalHarmonics[i - 1][i - 1] -
           (i - 1) * _sphericalHarmonics[i - 2][i - 2]) / (i);
}

void Harmonics::addComplex(real x, real y)
{
   std::complex<real> ephi1m(sqrt(2.0));
   std::complex<real> mult(x, y);

   for(size_t m = 1; m < _n; m++)
   {
      ephi1m *= mult;

      for(size_t l = m; l < _n; l++)
      {
         size_t center = _sphericalHarmonics[l].size() / 2;

         _sphericalHarmonics[l][center + m] *= ephi1m.real();
         _sphericalHarmonics[l][center - m] *= ephi1m.imag();
      }
   }
}

std::vector<std::vector<real>> Harmonics::calcSolidHarmonics(size_t n,
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

   for(size_t l = 0; l < n; l++)
   {
      size_t center = solidlHarmonics[l].size() / 2;

      for(size_t m = 0; m <= l; m++)
      {
         solidlHarmonics[l][center + m] *= curr * (isRegular ? 1.0 / _factorials[l + m] : _factorials[l - m]);
      }

      for(size_t m = 1; m <= l; m++)
      {
         solidlHarmonics[l][center - m] *= curr * (isRegular ? 1.0 / _factorials[l + m] : _factorials[l - m]);
      }

      curr *= -mult;
   }

   return solidlHarmonics;
}

real Harmonics::getHarmonic(int l, int m, const std::vector<std::vector<real>>& harmonics)
{
   return harmonics[l][harmonics[l].size() / 2 + m];
}
