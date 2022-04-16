#include <stdio.h>
#include <cstdlib>
#include <chrono>

#include "cuda_helper.hpp"

namespace cuda
{
   void tryKernelLaunch()
   {
      cudaError_t cudaStatus = cudaGetLastError();

      if(cudaStatus != cudaError::cudaSuccess)
      {
         printf("Kernel launch failed!\n");

         std::exit(0);
      }
   }

   void tryKernelSynchronize()
   {
      cudaError_t cudaStatus = cudaDeviceSynchronize();

      if(cudaStatus != cudaError::cudaSuccess)
      {
         printf("Kernel synchrinisation failed!\n");

         std::exit(0);
      }
   }
}