#pragma once
#include "harmonics.hpp"

class MultipoleTranslator
{
public:
   static ComplexHarmonicSeries translateMultipole(
      const ComplexHarmonicSeries& a,
      const ComplexHarmonicSeries& b);

   static RealHarmonicSeries translateMultipole(
      const RealHarmonicSeries& a,
      const RealHarmonicSeries& b);

   static HarmonicSeries<Vector3> translateMultipoleWithComplex(
      const HarmonicSeries<Vector3>& expansion,
      const Vector3& translation);

   static HarmonicSeries<Vector3> translateMultipoleWithReal(
      const HarmonicSeries<Vector3>& expansion,
      const Vector3& translation);

   static real multipoleTranslationFactor(int m, int mu);
};