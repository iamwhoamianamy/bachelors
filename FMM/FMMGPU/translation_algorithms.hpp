#pragma once
#include <vector>
#include "real.hpp"
#include "math.hpp"
#include "harmonics.hpp"

namespace kernels
{
   std::vector<real> addVectors(const std::vector<real>& a,
                                const std::vector<real>& b);

   void translateAllCPU(
      Vector3* result,
      const real* a,
      const Vector3* b,
      size_t harmonicCount,
      size_t harmonicOrder);

   void translateAllGPU(
      Vector3* result,
      const real* a,
      const Vector3* b,
      size_t harmonicCount,
      size_t harmonicOrder);

   void translateAllCPUMatrix(
      std::vector<Complex>& result,
      const std::vector<Complex>& a,
      const std::vector<Complex>& b,
      size_t harmonicCount,
      size_t harmonicOrder);

   void translateAllGPUMatrix(
      Complex* result,
      const Complex* a,
      const Complex* b,
      size_t harmonicCount,
      size_t harmonicOrder);
}
