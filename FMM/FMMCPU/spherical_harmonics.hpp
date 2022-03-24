#pragma once
#include <vector>
#include "real.hpp"
#include "vector3.hpp"

class SphericalHarmonics
{
private:
   int _n;
   std::vector<size_t> _factorials;
   std::vector<std::vector<real>> _sphericalHarmonics;
public:
   const int maxFactorialNum = 20;
   
   SphericalHarmonics(int n, const Vector3& vec);
   real getharmonic(int l, int m) const;

private:
   void calcFactorials();
   void calcSphericalHarmonics(const Vector3& vec);
   void initHarmonicArrays();
   void fillWithLegendrePolynomials(real z);
   void fillWithLegendrePolynomialDerivatives(real z);
   void mirrorLegendrePolynomialDerivatives(real z);
   real calcLegendrePolynomial(int i, real z);
   void addComplex(real x, real y);
};

