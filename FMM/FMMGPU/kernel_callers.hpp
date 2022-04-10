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
   void translateAllGPU(Vector3* result,
                        const real* a,
                        const Vector3* b,
                        size_t count, size_t order);
}
