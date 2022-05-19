#include "multipole_translator.hpp"

ComplexHarmonicSeries MultipoleTranslator::translateMultipole(
   const ComplexHarmonicSeries& a,
   const ComplexHarmonicSeries& b)
{
   ComplexHarmonicSeries res(a.order());

   for(int l = 0; l <= a.order(); l++)
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
                  res.getHarmonic(l, m) += a.getHarmonic(lambda, mu) *
                     b.getHarmonic(dl, dm) * multipoleTranslationFactor(m, mu);
               }
            }
         }
      }
   }

   return res;
}

RealHarmonicSeries MultipoleTranslator::translateMultipole(
   const RealHarmonicSeries& a,
   const RealHarmonicSeries& b)
{
   RealHarmonicSeries res(b.order());

   for(int l = 0; l <= b.order(); l++)
   {
      real zeroRes = 0;

      for(int lambda = 0; lambda <= l; lambda++)
      {
         int dl = l - lambda;

         for(int mu = -lambda; mu <= lambda; mu++)
         {
            if(-dl <= mu && mu <= +dl)
            {
               zeroRes += b.getHarmonic(lambda, mu) * a.getHarmonic(dl, mu) *
                  multipoleTranslationFactor(0, mu);
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

               real RR = b.getReal(lambda, mu);
               real IR = b.getImag(lambda, mu);

               real RM = 0;
               real IM = 0;

               real RnM = 0;
               real InM = 0;

               if(-dl <= dm && dm <= dl)
               {
                  RM = a.getReal(dl, dm);
                  IM = a.getImag(dl, dm);

                  realRes += (RR * RM - IR * IM) * multipoleTranslationFactor(m, mu);
                  imagRes += (RR * IM + IR * RM) * multipoleTranslationFactor(m, mu);
               }

               if(-dl <= dnm && dnm <= dl)
               {
                  RnM = a.getReal(dl, dnm);
                  InM = a.getImag(dl, dnm);

                  realRes += (RR * RnM - IR * InM) * multipoleTranslationFactor(-m, mu);
                  imagRes -= (RR * InM + IR * RnM) * multipoleTranslationFactor(-m, mu);
               }
            }
         }

         res.getHarmonic(l, m) = realRes * math::R_SQRT_2;
         res.getHarmonic(l, -m) = imagRes * math::R_SQRT_2;
      }
   }

   return res;
}

HarmonicSeries<Vector3> MultipoleTranslator::translateMultipoleWithComplex(
   const HarmonicSeries<Vector3>& expansion,
   const Vector3& translation)
{
   auto regular = Harmonics::realToComplex(Harmonics::calcRegularSolidHarmonics(expansion.order(), translation));

   auto xComponent = Harmonics::realToComplex(Harmonics::separateCoord(expansion, 0));
   auto yComponent = Harmonics::realToComplex(Harmonics::separateCoord(expansion, 1));
   auto zComponent = Harmonics::realToComplex(Harmonics::separateCoord(expansion, 2));

   return Harmonics::createFormXYZ(
      Harmonics::complexToReal(translateMultipole(regular, xComponent)),
      Harmonics::complexToReal(translateMultipole(regular, yComponent)),
      Harmonics::complexToReal(translateMultipole(regular, zComponent)));
}

HarmonicSeries<Vector3> MultipoleTranslator::translateMultipoleWithReal(
   const HarmonicSeries<Vector3>& expansion,
   const Vector3& translation)
{
   auto regular = Harmonics::calcRegularSolidHarmonics(expansion.order(), translation);

   auto xComponent = Harmonics::separateCoord(expansion, 0);
   auto yComponent = Harmonics::separateCoord(expansion, 1);
   auto zComponent = Harmonics::separateCoord(expansion, 2);

   return Harmonics::createFormXYZ(
      translateMultipole(regular, xComponent),
      translateMultipole(regular, yComponent),
      translateMultipole(regular, zComponent));
}

real MultipoleTranslator::multipoleTranslationFactor(int m, int mu)
{
   return pow(-1, -0.5 * (abs(m) - abs(mu) - abs(m - mu)));
}

ComplexHarmonicSeries MultipoleTranslator::multipoleToLocal(
   const ComplexHarmonicSeries& a,
   const ComplexHarmonicSeries& b)
{
   ComplexHarmonicSeries res(a.order());

   for(int l = 0; l <= a.order(); l++)
   {
      for(int lambda = 0; lambda <= a.order(); lambda++)
      {
         for(int m = -l; m <= l; m++)
         {
            for(int mu = -lambda; mu <= lambda; mu++)
            {
               int dl = l + lambda;
               int dm = m - mu;

               if(dm >= -dl && dm <= +dl && dl <= a.order())
               {
                  res.getHarmonic(l, m) += a.getHarmonic(lambda, mu) *
                     b.getHarmonic(dl, dm) * multipoleToLocalTranslationFactor(m, mu, lambda);
               }
            }
         }
      }
   }

   return res;
}

HarmonicSeries<Vector3> MultipoleTranslator::multipoleToLocalWithComplex(
   const HarmonicSeries<Vector3>& expansion,
   const Vector3& translation)
{
   auto irregular = Harmonics::realToComplex(
      Harmonics::calcIrregularSolidHarmonics(expansion.order(), translation));

   auto xComponent = Harmonics::realToComplex(Harmonics::separateCoord(expansion, 0));
   auto yComponent = Harmonics::realToComplex(Harmonics::separateCoord(expansion, 1));
   auto zComponent = Harmonics::realToComplex(Harmonics::separateCoord(expansion, 2));

   return Harmonics::createFormXYZ(
      Harmonics::complexToReal(multipoleToLocal(irregular, xComponent)),
      Harmonics::complexToReal(multipoleToLocal(irregular, yComponent)),
      Harmonics::complexToReal(multipoleToLocal(irregular, zComponent)));

   //return Harmonics::createFormXYZ(
   //   Harmonics::complexToReal(multipoleToLocal(xComponent, irregular)),
   //   Harmonics::complexToReal(multipoleToLocal(yComponent, irregular)),
   //   Harmonics::complexToReal(multipoleToLocal(zComponent, irregular)));
}

real MultipoleTranslator::multipoleToLocalTranslationFactor(int m, int mu, int lambda)
{
   return pow(-1, -0.5 * (abs(m - mu) - abs(m) - abs(mu))) * 
      pow(-1, abs(lambda));
}

ComplexHarmonicSeries MultipoleTranslator::translateLocal(
   const ComplexHarmonicSeries& a,
   const ComplexHarmonicSeries& b)
{
   ComplexHarmonicSeries res(a.order());

   for(int l = 0; l <= a.order(); l++)
   {
      for(int lambda = 0; lambda <= a.order(); lambda++)
      {
         for(int m = -l; m <= l; m++)
         {
            for(int mu = -lambda; mu <= lambda; mu++)
            {
               int dl = lambda - l;
               int dm = m - mu;

               if(dm >= -dl && dm <= +dl && dl <= a.order())
               {
                  res.getHarmonic(l, m) += a.getHarmonic(lambda, mu) *
                     b.getHarmonic(dl, dm) * localTranslationFactor(m, mu, lambda, l);
               }
            }
         }
      }
   }

   return res;
}

HarmonicSeries<Vector3> MultipoleTranslator::translateLocalWithComplex(
   const HarmonicSeries<Vector3>& expansion,
   const Vector3& translation)
{
   auto regular = Harmonics::realToComplex(
      Harmonics::calcRegularSolidHarmonics(expansion.order(), translation));

   auto xComponent = Harmonics::realToComplex(Harmonics::separateCoord(expansion, 0));
   auto yComponent = Harmonics::realToComplex(Harmonics::separateCoord(expansion, 1));
   auto zComponent = Harmonics::realToComplex(Harmonics::separateCoord(expansion, 2));

   return Harmonics::createFormXYZ(
      Harmonics::complexToReal(translateLocal(regular, xComponent)),
      Harmonics::complexToReal(translateLocal(regular, yComponent)),
      Harmonics::complexToReal(translateLocal(regular, zComponent)));

   //return Harmonics::createFormXYZ(
   //   Harmonics::complexToReal(translateLocal(xComponent, regular)),
   //   Harmonics::complexToReal(translateLocal(yComponent, regular)),
   //   Harmonics::complexToReal(translateLocal(zComponent, regular)));
}

real MultipoleTranslator::localTranslationFactor(int m, int mu, int lambda, int l)
{
   return pow(-1, -0.5 * (abs(mu) - abs(m - mu) - abs(m))) *
      pow(-1, abs(lambda) + abs(l));
}