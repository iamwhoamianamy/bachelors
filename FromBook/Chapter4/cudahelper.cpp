#pragma once
#include "cudahelper.h"
#include "cuda_runtime.h"
#include <cstdlib>
#include <stdio.h>

cudaError_t HANDLE_ERROR(cudaError_t funct_return)
{
   switch(funct_return)
   {
      case cudaError::cudaSuccess:
      {
         return funct_return;
      }
      default:
      {
         printf("ERROR!\n");
         std::exit(0);
      }
   }
}

cudaError_t TryCudaMalloc(void** dev_ptr, const size_t size)
{
   cudaError_t cuda_status = cudaMalloc(dev_ptr, size);

   if(cuda_status != cudaError::cudaSuccess)
   {
      printf("ERROR!\n");

      free(*dev_ptr);

      std::exit(0);
   }
   else
      return cuda_status;
}

