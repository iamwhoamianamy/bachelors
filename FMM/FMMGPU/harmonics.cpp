#include <thrust/complex.h>
#include <limits>
#include "harmonics.hpp"
#include "math.hpp"
#include <iostream>

Factorials Harmonics::_factorials;

Harmonics::Harmonics(size_t order, const Vector3& point) :
   _order(order)
{
   calcSphericalHarmonics(point);
}

const RealHarmonicSeries& Harmonics::sphericalHarmonics() const
{
   return _sphericalHarmonics;
}

RealHarmonicSeries Harmonics::calcSolidHarmonics(size_t order,
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

   auto solidlHarmonics = Harmonics(order, point).sphericalHarmonics();

   real mult = isRegular ? r : 1 / r;
   real curr = isRegular ? 1 : mult;

   for(int l = 0; l <= order; l++)
   {
      for(int m = 0; m <= l; m++)
      {
         solidlHarmonics.getHarmonic(l, -m) *= curr *
            (isRegular ? 1.0 / _factorials[l + m] : _factorials[l - m]);
      }

      for(int m = 1; m <= l; m++)
      {
         solidlHarmonics.getHarmonic(l, m) *= curr *
            (isRegular ? 1.0 / _factorials[l + m] : _factorials[l - m]);
      }

      curr *= -mult;
   }

   return solidlHarmonics;
}

RealHarmonicSeries Harmonics::calcRegularSolidHarmonics(size_t order, Vector3 point)
{
   return calcSolidHarmonics(order, point, true);
}

RealHarmonicSeries Harmonics::calcIrregularSolidHarmonics(size_t order, Vector3 point)
{
   return calcSolidHarmonics(order, point, false);
}

void Harmonics::calcSphericalHarmonics(const Vector3& point)
{
   _sphericalHarmonics = RealHarmonicSeries(_order);
   fillWithLegendrePolynomials(point.z);
   fillWithLegendrePolynomialDerivatives(point.z);
   mirrorLegendrePolynomialDerivatives(point.z);
   addComplex(point.x, point.y);
}

void Harmonics::fillWithLegendrePolynomials(real z)
{
   _sphericalHarmonics.getHarmonic(0, 0) = 1;
   _sphericalHarmonics.getHarmonic(1, 0) = z;

   for(size_t l = 2; l <= _order; l++)
   {
      _sphericalHarmonics.getHarmonic(l, 0) = calcLegendrePolynomial(l, z);
   }
}

void Harmonics::fillWithLegendrePolynomialDerivatives(real z)
{
   _sphericalHarmonics.getHarmonic(0, 0) = 1;
   _sphericalHarmonics.getHarmonic(1, -1) = 1;
   _sphericalHarmonics.getHarmonic(1, 1) = 1;

   for(size_t l = 2; l <= _order; l++)
   {
      for(size_t m = 1; m < l; m++)
      {
         _sphericalHarmonics.getHarmonic(l, m) = z * _sphericalHarmonics.getHarmonic(l - 1, m) +
            (l + m - 1) * _sphericalHarmonics.getHarmonic(l - 1, m - 1);
      }

      _sphericalHarmonics.getHarmonic(l, l) = (2 * l - 1) *
         _sphericalHarmonics.getHarmonic(l - 1, l - 1);
   }

}

void Harmonics::mirrorLegendrePolynomialDerivatives(real z)
{
   for(int l = 2; l <= _order; l++)
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
   Complex ephi1m = makeComplex(sqrt(2.0), 0);
   Complex mult = makeComplex(x, y);

   for(int m = 1; m < _order; m++)
   {
      ephi1m *= mult;

      for(int l = m; l <= _order; l++)
      {
         _sphericalHarmonics.getHarmonic(l, m) *= getReal(ephi1m);
         _sphericalHarmonics.getHarmonic(l, -m) *= getImag(ephi1m);
      }
   }
}

ComplexMatrix Harmonics::calcRealToComplexMatrix2D(size_t order)
{
   size_t harmonicLen = (order + 1) * (order + 1);
   ComplexMatrix res(harmonicLen, std::vector<Complex>(harmonicLen));

   for(int l = 0; l <= order; l++)
   {
      for(int m = -l; m <= l; m++)
      {
         real realPart = HarmonicSeries<real>::getRealFactor(l, m);
         real imagPart = HarmonicSeries<real>::getImagFactor(l, m);

         int realM = abs(m);
         int imagM = -abs(m);

         size_t currentIndex = HarmonicSeries<Complex>::lmToIndex(l, m);
         size_t newRealIndex = HarmonicSeries<Complex>::lmToIndex(l, realM);
         size_t newImagIndex = HarmonicSeries<Complex>::lmToIndex(l, imagM);

         res[newRealIndex][currentIndex] += makeComplex(realPart, 0);
         res[newImagIndex][currentIndex] += makeComplex(0, imagPart);
      }
   }

   return res;
}

ComplexMatrix Harmonics::calcComplexToRealMatrix2D(size_t order)
{
   size_t harmonicLen = (order + 1) * (order + 1);
   ComplexMatrix res(harmonicLen, std::vector<Complex>(harmonicLen));

   for(int l = 0; l <= order; l++)
   {
      res[l * l + l][l * l + l] = makeComplex(1, 0);

      for(int m = 1; m <= l; m++)
      {
         res[l * l + l + m][l * l + l + m] = makeComplex(math::R_SQRT_2, 0);
         res[l * l + l - m][l * l + l + m] = makeComplex(math::R_SQRT_2, 0);

         res[l * l + l + m][l * l + l - m] = makeComplex(0, -math::R_SQRT_2);
         res[l * l + l - m][l * l + l - m] = makeComplex(0, math::R_SQRT_2);
      }
   }

   return res;
}

std::vector<Complex> Harmonics::calcRealToComplexTransitionMatrix1D(size_t order)
{
   size_t harmonicLen = (order + 1) * (order + 1);
   std::vector<Complex> res(harmonicLen * harmonicLen);

   for(int l = 0; l <= order; l++)
   {
      for(int m = -l; m <= l; m++)
      {
         real realPart = HarmonicSeries<real>::getRealFactor(l, m);
         real imagPart = HarmonicSeries<real>::getImagFactor(l, m);

         int realM = abs(m);
         int imagM = -abs(m);

         size_t currentIndex = HarmonicSeries<Complex>::lmToIndex(l, m);
         size_t newRealIndex = HarmonicSeries<Complex>::lmToIndex(l, realM);
         size_t newImagIndex = HarmonicSeries<Complex>::lmToIndex(l, imagM);

         res[newRealIndex * harmonicLen + currentIndex] +=
            makeComplex(realPart, 0);
         res[newImagIndex * harmonicLen + currentIndex] +=
            makeComplex(0, imagPart);
      }
   }

   return res;
}

std::vector<Complex> Harmonics::calcComplexToRealTransitionMatrix1D(size_t order)
{
   size_t harmonicLen = (order + 1) * (order + 1);
   std::vector<Complex> res(harmonicLen * harmonicLen);

   for(int l = 0; l <= order; l++)
   {
      res[(l * l + l) * harmonicLen + (l * l + l)] = makeComplex(1, 0);

      for(int m = 1; m <= l; m++)
      {
         res[(l * l + l + m) * harmonicLen + (l * l + l + m)] =
            makeComplex(math::R_SQRT_2, 0);
         res[(l * l + l - m) * harmonicLen + (l * l + l + m)] =
            makeComplex(math::R_SQRT_2, 0);

         res[(l * l + l + m) * harmonicLen + (l * l + l - m)] =
            makeComplex(0, -math::R_SQRT_2);
         res[(l * l + l - m) * harmonicLen + (l * l + l - m)] =
            makeComplex(0, math::R_SQRT_2);
      }
   }

   return res;
}

real Harmonics::getFactorial(size_t n)
{
   return _factorials[n];
}

RealHarmonicSeries Harmonics::separateCoord(
   const HarmonicSeries<Vector3>& harmonics,
   size_t i)
{
   RealHarmonicSeries res(harmonics.order());

   for(size_t h = 0; h < harmonics.elemCount(); h++)
   {
      res.getHarmonic(h) = harmonics.getHarmonic(h)[i];
   }

   return res;
}

HarmonicSeries<Vector3> Harmonics::createFormXYZ(const RealHarmonicSeries& xs,
                                                 const RealHarmonicSeries& ys,
                                                 const RealHarmonicSeries& zs)
{
   HarmonicSeries<Vector3> res(xs.order());

   for(size_t i = 0; i < xs.elemCount(); i++)
   {
      res.getHarmonic(i) = Vector3(xs.getHarmonic(i), 
                                   ys.getHarmonic(i), 
                                   zs.getHarmonic(i));
   }

   return res;
}

ComplexHarmonicSeries Harmonics::realToComplex(
   const RealHarmonicSeries& harmonics)
{
   ComplexHarmonicSeries res(harmonics.order());
   
   real c = math::R_SQRT_2;
   for(int l = 0; l <= harmonics.order(); l++)
   {
      res.getHarmonic(l, 0) = makeComplex(harmonics.getHarmonic(l, 0), 0);
      for(int m = 1; m <= l; m++)
      {
         res.setHarmonic(l, m, makeComplex(
            harmonics.getHarmonic(l, m) * c,
            harmonics.getHarmonic(l, -m) * c));
         res.setHarmonic(l, -m, makeComplex(
            harmonics.getHarmonic(l, m) * c,
            -harmonics.getHarmonic(l, -m) * c));
      }
   }

   return res;
}

RealHarmonicSeries Harmonics::complexToReal(const ComplexHarmonicSeries& harmonics)
{
   RealHarmonicSeries res(harmonics.order());

   real c = 1.0 / sqrt(2);
   for(int l = 0; l <= harmonics.order(); l++)
   {
      res.getHarmonic(l, 0) = getReal(harmonics.getHarmonic(l, 0));
      for(int m = 1; m <= l; m++)
      {
         res.getHarmonic(l, m) = c * getReal(harmonics.getHarmonic(l, m) + harmonics.getHarmonic(l, -m));
         res.getHarmonic(l, -m) = c * getImag(harmonics.getHarmonic(l, m) - harmonics.getHarmonic(l, -m));
      }
   }

   return res;
}

