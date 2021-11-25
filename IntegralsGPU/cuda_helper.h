#pragma once
#include "real.h"
#include "cuda_runtime.h"
//#include <windows.h>

namespace cuda_utilities
{
   typedef void handler;
#define FOREGROUND_WHITE 7

   //void printg(const char* str)
   //{
   //   SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_GREEN);
   //   printf(str);
   //   SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_WHITE);
   //}

   //void printr(const char* str)
   //{
   //   SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_RED);
   //   printf(str);
   //   SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_WHITE);
   //}

   //cudaError_t handleError(cudaError_t functReturn)
   //{
   //   switch(functReturn)
   //   {
   //      case cudaError::cudaSuccess:
   //      {
   //         return functReturn;
   //      }
   //      default:
   //      {
   //         printr("ERROR!\n");
   //         std::exit(0);
   //      }
   //   }
   //}

   //void tryCudaMalloc(void** devPtr, const size_t size)
   //{
   //   cudaError_t cudaStatus = cudaMalloc(devPtr, size);

   //   if(cudaStatus != cudaError::cudaSuccess)
   //   {
   //      printr("Error on allocating memory!\n");

   //      free(*devPtr);

   //      std::exit(0);
   //   }
   //}

   //void tryCudaLastError()
   //{
   //   cudaError_t cudaStatus = cudaGetLastError();

   //   if(cudaStatus != cudaError::cudaSuccess)
   //   {
   //      printr("An error occured!\n");

   //      std::exit(0);
   //   }
   //}

   //void tryCudaReset()
   //{
   //   cudaError_t cudaStatus = cudaDeviceReset();

   //   if(cudaStatus == cudaError::cudaSuccess)
   //   {
   //      printg("Cuda device reset success!\n");
   //   }
   //   else
   //   {
   //      printr("Cuda device reset failed!\n");

   //      std::exit(0);
   //   }
   //}

   void tryKernelLaunch();
   void tryKernelSynchronize();

   // Gets next closest integer to number that is devisible on devider
   size_t nextDevisible(const size_t number, const size_t devidor);
}