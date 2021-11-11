#include <stdio.h>
#include <cstdlib>
#include "cuda_helper.h"

void cuda_utilities::tryKernelLaunch()
{
   cudaError_t cudaStatus = cudaGetLastError();

   if(cudaStatus != cudaError::cudaSuccess)
   {
      printf("Kernel launch failed!\n");

      std::exit(0);
   }
}

void cuda_utilities::tryKernelSynchronize()
{
   cudaError_t cudaStatus = cudaDeviceSynchronize();

   if(cudaStatus != cudaError::cudaSuccess)
   {
      printf("Kernel synchrinisation failed!\n");

      std::exit(0);
   }
}