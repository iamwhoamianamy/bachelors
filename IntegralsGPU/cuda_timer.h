#pragma once
#include"real.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

namespace cuda_utilities
{
   class CudaTimer
   {
   private:
      cudaEvent_t start, stop;

   public:
      CudaTimer();
      ~CudaTimer();
      void Start();
      real Ellapsed();
   };
}