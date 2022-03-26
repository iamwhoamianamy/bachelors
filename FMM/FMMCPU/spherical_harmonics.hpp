#pragma once
#include <vector>
#include "real.hpp"
#include "vector3.hpp"
#include "factorials.hpp"

class SphericalHarmonics
{
private:
   int _n;
   std::vector<std::vector<real>> _sphericalHarmonics;
   static Factorials _factorials;
public:   
   SphericalHarmonics(int n, const Vector3& point);
   real getHarmonic(int l, int m) const;

   const std::vector<std::vector<real>>& sphericalHarmonics() const;
   static std::vector<std::vector<real>> calcSolidHarmonics(size_t n, 
                                                            Vector3 point,
                                                            bool isRegular);

private:
   void calcSphericalHarmonics(const Vector3& point);
   void initHarmonicArrays();
   void fillWithLegendrePolynomials(real z);
   void fillWithLegendrePolynomialDerivatives(real z);
   void mirrorLegendrePolynomialDerivatives(real z);
   real calcLegendrePolynomial(int i, real z);
   void addComplex(real x, real y);
};

