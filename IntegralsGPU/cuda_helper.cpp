#include <stdio.h>
#include <cstdlib>
#include <chrono>

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

size_t cuda_utilities::nextDevisible(const size_t number, const size_t devidor)
{
   if(devidor == 0)
      return number;
   else
      return (number + devidor - 1) / devidor * devidor;
}
