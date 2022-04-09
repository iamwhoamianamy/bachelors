#pragma once
#include "cuda_runtime.h"

namespace cuda
{
   void tryKernelLaunch();
   void tryKernelSynchronize();

   // Gets next closest integer to number that is devisible on devider
   size_t nextDevisible(const size_t number, const size_t devidor);
}