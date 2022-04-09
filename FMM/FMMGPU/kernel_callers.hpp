#pragma once
#include <vector>
#include "real.hpp"
#include "math.hpp"
#include "harmonics.hpp"

namespace kernels
{
   std::vector<real> addVectors(const std::vector<real>& a,
                                const std::vector<real>& b);

   void translateAllCPU(Vector3* result,
                        const real* a,
                        const Vector3* b,
                        size_t count, size_t order);

   template <typename T>
   T* getHarmonic(T* harmonics,
                 int harmonicBegin,
                 int l, int m)
   {
      return &harmonics[harmonicBegin + l * l + l + m];
   }

   template <typename T>
   T getReal(const T* harmonics,
             int harmonicBegin,
             int l, int m)
   {
      return *getHarmonic(harmonics, harmonicBegin, l, abs(m)) *
         (math::R_SQRT_2 * (m != 0) + (m == 0));
   }

   template <typename T>
   T getImag(const T* harmonics,
             int harmonicBegin,
             int l, int m)
   {
      return *getHarmonic(harmonics, harmonicBegin, l, -abs(m)) * 
         math::R_SQRT_2 * (m != 0) * math::sign(m);
   }
}
