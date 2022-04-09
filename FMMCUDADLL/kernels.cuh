#include "cuda_runtime.h"
#include "device_launch_parameters.h"

namespace kernels
{
   template<class T>
   __global__ void addingKernel(T* c, const T* a, const T* b)
   {
      unsigned int i = threadIdx.x;
      c[i] = a[i] + b[i];
   }
}