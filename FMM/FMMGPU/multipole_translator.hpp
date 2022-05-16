#pragma once
#include "harmonics.hpp"

class MultipoleTranslator
{
public:

#pragma region Multipole to multipole algorithms

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

#pragma endregion

#pragma region Multipole to local algorithms

   static ComplexHarmonicSeries multipoleToLocal(
      const ComplexHarmonicSeries& a,
      const ComplexHarmonicSeries& b);

   static HarmonicSeries<Vector3> multipoleToLocalWithComplex(
      const HarmonicSeries<Vector3>& expansion,
      const Vector3& translation);

   static real multipoleToLocalTranslationFactor(int m, int mu, int lambda);

#pragma endregion
};