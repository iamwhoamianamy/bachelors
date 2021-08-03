#pragma once
#include"cuda_runtime.h"
#include<cstdlib>
#include <stdio.h>

int HANDLE_ERROR(cudaError_t funct_return)
{
   switch(funct_return)
   {
      case cudaError::cudaSuccess:
      {
         return 0;
      }
      default:
      {
         printf("ERROR!\n");
         std::exit(0);
      }
   }
}