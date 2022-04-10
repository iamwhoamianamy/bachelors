#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "harmonics.hpp"

namespace kernels
{
   template<class T>
   __global__ void addingKernel(T* c, const T* a, const T* b)
   {
      uint i = threadIdx.x;
      c[i] = a[i] + b[i];
   }

   __global__ void translateAllGPUKernelSimple(
      Vector3* result, const real* a, const Vector3* b,
      size_t count, size_t order);

   __global__ void translateAllGPUKernelBlockForHarmonic(
      Vector3* result, const real* a, const Vector3* b,
      size_t count, size_t order);

   __all__ size_t lmToIndex(int harmonicBegin,
                    int l, int m);

   template <typename T>
   __all__ T& getHarmonic(T* harmonics,
                  int harmonicBegin,
                  int l, int m)
   {
      return harmonics[harmonicBegin + l * l + l + m];
   }

   template <typename T>
   __all__ T getReal(const T* harmonics,
             int harmonicBegin,
             int l, int m)
   {
      return getHarmonic(harmonics, harmonicBegin, l, abs(m)) *
         (math::R_SQRT_2 * (m != 0) + (m == 0));
   }

   template <typename T>
   __all__ T getImag(const T* harmonics,
             int harmonicBegin,
             int l, int m)
   {
     return getHarmonic(harmonics, harmonicBegin, l, -abs(m)) *
         math::R_SQRT_2 * (m != 0) * math::sign(m);
   }

   __all__ real strangeFactor(int m, int mu);
}