#include <complex>
#include <limits>
#include "harmonics.hpp"
#include "math.hpp"

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
         _sphericalHarmonics.getHarmonic(l, m) *= ephi1m.real();
         _sphericalHarmonics.getHarmonic(l, -m) *= ephi1m.imag();
      }
   }
}

real Harmonics::strangeFactor(int m, int mu)
{
   return pow(-1, -0.5 * (abs(m) - abs(mu) - abs(m - mu)));
}

HarmonicSeries<std::complex<real>> Harmonics::translate(
   int n,
   const HarmonicSeries<std::complex<real>>& regular,
   const HarmonicSeries<std::complex<real>>& child)
{
   HarmonicSeries<std::complex<real>> res(n);

   for(int l = 0; l < n; l++)
   {
      for(int m = -l; m <= l; m++)
      {
         for(int lambda = 0; lambda <= l; lambda++)
         {
            int dl = l - lambda;

            for(int mu = -lambda; mu <= lambda; mu++)
            {
               int dm = m - mu;

               if(dm >= -dl && dm <= +dl)
               {
                  res.getHarmonic(l, m) += regular.getHarmonic(lambda, mu) *
                     child.getHarmonic(dl, dm) * strangeFactor(m, mu);
               }
            }
         }
      }
   }

   return res;
}

HarmonicSeries<real> Harmonics::translate(
   int n,
   const HarmonicSeries<real>& regular,
   const HarmonicSeries<real>& child)
{
   HarmonicSeries<real> res(n);

   for(int l = 0; l < n; l++)
   {
      real zeroRes = 0;

      for(int lambda = 0; lambda <= l; lambda++)
      {
         int dl = l - lambda;

         for(int mu = -lambda; mu <= lambda; mu++)
         {
            if(-dl <= mu && mu <= +dl)
            {
               zeroRes += regular.getHarmonic(lambda, mu) * child.getHarmonic(dl, mu) *
                           strangeFactor(0, mu);
            }
         }
      }

      res.getHarmonic(l, 0) = zeroRes;

      for(int m = 1; m <= l; m++)
      {
         real realRes = 0;
         real imagRes = 0;

         for(int lambda = 0; lambda <= l; lambda++)
         {
            int dl = l - lambda;

            for(int mu = -lambda; mu <= lambda; mu++)
            {
               int dm = m - mu;
               int dnm = -m - mu;

               real RR = regular.getReal(lambda, mu);
               real IR = regular.getImag(lambda, mu);

               real RM = 0;
               real IM = 0;

               real RnM = 0;
               real InM = 0;

               if(-dl <= dm && dm <= dl)
               {
                  RM = child.getReal(dl, dm);
                  IM = child.getImag(dl, dm);

                  realRes += (RR * RM - IR * IM) * strangeFactor(m, mu);
                  imagRes += (RR * IM + IR * RM) * strangeFactor(m, mu);
               }

               if(-dl <= dnm && dnm <= dl)
               {
                  RnM = child.getReal(dl, dnm);
                  InM = child.getImag(dl, dnm);

                  realRes += (RR * RnM - IR * InM) * strangeFactor(-m, mu);
                  imagRes -= (RR * InM + IR * RnM) * strangeFactor(-m, mu);
               }
            }
         }

         res.getHarmonic(l, m) = realRes * math::R_SQRT_2;
         res.getHarmonic(l, -m) = imagRes * math::R_SQRT_2;
      }
   }

   return res;
}

real Harmonics::getRealChildContribution(int l, int m,
                                         const HarmonicSeries<real>& regular,
                                         const HarmonicSeries<real>& child)
{
   real res = 0;

   for(int lambda = 0; lambda <= l; lambda++)
   {
      int dl = l - lambda;

      for(int mu = -lambda; mu <= lambda; mu++)
      {
         int dm = m - mu;

         if(dm >= -dl && dm <= +dl)
         {
            res += (regular.getReal(lambda, mu) * child.getReal(dl, dm) -
                    regular.getImag(lambda, mu) * child.getImag(dl, dm)) * 
                    strangeFactor(m, mu);
         }
      }
   }

   return res;
}

real Harmonics::getImagChildContribution(int l, int m,
                                         const HarmonicSeries<real>& regular,
                                         const HarmonicSeries<real>& child)
{
   real res = 0;

   for(int lambda = 0; lambda <= l; lambda++)
   {
      int dl = l - lambda;

      for(int mu = -lambda; mu <= lambda; mu++)
      {
         int dm = m - mu;

         if(dm >= -dl && dm <= +dl)
         {
            res += (regular.getReal(lambda, mu) * child.getImag(dl, dm) +
                    regular.getImag(lambda, mu) * child.getReal(dl, dm)) * 
                    strangeFactor(m, mu);
         }
      }
   }

   return res;
}

real Harmonics::getFactorial(size_t n)
{
   return _factorials[n];
}

HarmonicSeries<real> Harmonics::separateX(const HarmonicSeries<Vector3>& harmonics)
{
   HarmonicSeries<real> res(harmonics.size());

   for(int l = 0; l < harmonics.size(); l++)
   {
      for(int m = -l; m <= l; m++)
      {
         res.getHarmonic(l, m) = harmonics.getHarmonic(l, m).x;
      }
   }

   return res;
}

HarmonicSeries<real> Harmonics::separateY(const HarmonicSeries<Vector3>& harmonics)
{
   HarmonicSeries<real> res(harmonics.size());

   for(int l = 0; l < harmonics.size(); l++)
   {
      for(int m = -l; m <= l; m++)
      {
         res.getHarmonic(l, m) = harmonics.getHarmonic(l, m).y;
      }
   }

   return res;
}

HarmonicSeries<real> Harmonics::separateZ(const HarmonicSeries<Vector3>& harmonics)
{
   HarmonicSeries<real> res(harmonics.size());

   for(int l = 0; l < harmonics.size(); l++)
   {
      for(int m = -l; m <= l; m++)
      {
         res.getHarmonic(l, m) = harmonics.getHarmonic(l, m).z;
      }
   }

   return res;
}

HarmonicSeries<Vector3> Harmonics::createFormXYZ(const HarmonicSeries<real>& xs,
                                                 const HarmonicSeries<real>& ys,
                                                 const HarmonicSeries<real>& zs)
{
   HarmonicSeries<Vector3> res(xs.size());

   for(int l = 0; l < xs.size(); l++)
   {
      for(int m = -l; m <= l; m++)
      {
         res.getHarmonic(l, m) = Vector3(xs.getHarmonic(l, m), ys.getHarmonic(l, m), zs.getHarmonic(l, m));
      }
   }

   return res;
}

HarmonicSeries<std::complex<real>> Harmonics::realToComplex(
   const HarmonicSeries<real>& harmonics)
{
   HarmonicSeries<std::complex<real>> res(harmonics.size());
   
   real c = 1.0 / sqrt(2);
   for(int l = 0; l < harmonics.size(); l++)
   {
      res.getHarmonic(l, 0) = std::complex<real>(harmonics.getHarmonic(l, 0), 0);
      for(int m = 1; m <= l; m++)
      {
         res.getHarmonic(l, m) =  c * std::complex<real>(harmonics.getHarmonic(l, m),
                                                         harmonics.getHarmonic(l, -m));
         res.getHarmonic(l, -m) = c * std::complex<real>(harmonics.getHarmonic(l, m),
                                                        -harmonics.getHarmonic(l, -m));
      }
   }

   return res;
}

HarmonicSeries<real> Harmonics::complexToReal(const HarmonicSeries<std::complex<real>>& harmonics)
{
   HarmonicSeries<real> res(harmonics.size());

   real c = 1.0 / sqrt(2);
   for(int l = 0; l < harmonics.size(); l++)
   {
      res.getHarmonic(l, 0) = harmonics.getHarmonic(l, 0).real();
      for(int m = 1; m <= l; m++)
      {
         res.getHarmonic(l, m) = c *  (harmonics.getHarmonic(l, m) + harmonics.getHarmonic(l, -m)).real();
         res.getHarmonic(l, -m) = c * (harmonics.getHarmonic(l, m) - harmonics.getHarmonic(l, -m)).imag();
      }
   }

   return res;
}

