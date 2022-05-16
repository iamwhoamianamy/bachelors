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
      Harmonics::complexToReal(
         MultipoleTranslator::translateMultipole(regular, xComponent)),
      Harmonics::complexToReal(
         MultipoleTranslator::translateMultipole(regular, yComponent)),
      Harmonics::complexToReal(
         MultipoleTranslator::translateMultipole(regular, zComponent)));
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
      MultipoleTranslator::translateMultipole(regular, xComponent),
      MultipoleTranslator::translateMultipole(regular, yComponent),
      MultipoleTranslator::translateMultipole(regular, zComponent));
}

real MultipoleTranslator::multipoleTranslationFactor(int m, int mu)
{
   return pow(-1, -0.5 * (abs(m) - abs(mu) - abs(m - mu)));
}