#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "device_functions.h"
#include "cudahelper.cuh"
#include <stdio.h>
#include <memory>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <chrono>

const int N = 256 * 1024 * 1024;            // Arrays size
const int threadsPerBlock = 256;    // Emount of participating threads per block
const int maxbBlocksPerGrid = 32;

// Emount of participating blocks
//const int blocksPerGrid = min(maxbBlocksPerGrid, (N + threadsPerBlock - 1) / threadsPerBlock);
const int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;

__global__ void dot1(double* a, double* b, double* c)
{
   __shared__ double cache[threadsPerBlock];
   int tid = threadIdx.x + blockIdx.x * blockDim.x;
   int cacheIndex = threadIdx.x;

   double temp = 0;

   while(tid < N)
   {
      temp += a[tid] * b[tid];
      tid += blockDim.x * gridDim.x;
   }

   cache[cacheIndex] = temp;

   __syncthreads();

   int i = threadsPerBlock / 2;

   while(i != 0)
   {
      if(cacheIndex < i)
         cache[cacheIndex] += cache[cacheIndex + i];
      __syncthreads();
      i /= 2;
   }

   //int numOfSums = threadsPerBlock / 2;

   //while(numOfSums > 0)
   //{
   //   if(cacheIndex < numOfSums)
   //   {
   //      int oldIndex = cacheIndex * 2;
   //      cache[cacheIndex] = cache[oldIndex] + cache[oldIndex + 1];
   //   }
   //   __syncthreads();
   //   numOfSums /= 2;
   //}

   if(cacheIndex == 0)
      c[blockIdx.x] = cache[0];
}

double sumOfSquares(int n)
{
   return (double)n * (n + 1.0) * (2.0 * n + 1.0) / 6.0;
}

double sumOfN(int n)
{
   return n / 2.0 * (n - 1);
}

double testDot(_Out_ float &time)
{
   double* a, *b, *part_c;
   double* dev_a, * dev_b, * dev_part_c;

   a = (double*)std::malloc(N * sizeof(double));
   b = (double*)std::malloc(N * sizeof(double));
   part_c = (double*)std::malloc(blocksPerGrid * sizeof(double));

   tryCudaMalloc((void**)&dev_a, N * sizeof(double));
   tryCudaMalloc((void**)&dev_b, N * sizeof(double));
   tryCudaMalloc((void**)&dev_part_c, blocksPerGrid * sizeof(double));

   for(size_t i = 0; i < N; i++)
   {
      //a[i] = i;
      //b[i] = i * 2;

      a[i] = 1;
      b[i] = i;
   }

   handleError(cudaMemcpy(dev_a, a, N * sizeof(double), cudaMemcpyHostToDevice));
   handleError(cudaMemcpy(dev_b, b, N * sizeof(double), cudaMemcpyHostToDevice));
   
   cudaEvent_t start, stop;
   handleError(cudaEventCreate(&start));
   handleError(cudaEventCreate(&stop));
   handleError(cudaEventRecord(start, 0));

   dot1 << <blocksPerGrid, threadsPerBlock >> > (dev_a, dev_b, dev_part_c);

   tryKernelLaunch();
   tryKernelSynchronize();

   handleError(cudaEventRecord(stop, 0));
   handleError(cudaEventSynchronize(stop));

   handleError(cudaEventElapsedTime(&time, start, stop));

   handleError(cudaEventDestroy(start));
   handleError(cudaEventDestroy(stop));

   handleError(cudaMemcpy(part_c, dev_part_c, blocksPerGrid * sizeof(double), cudaMemcpyDeviceToHost));

   cudaFree(dev_a);
   cudaFree(dev_b);
   cudaFree(dev_part_c);

   double dot_product = 0;

   for(size_t i = 0; i < blocksPerGrid; i++)
      dot_product += part_c[i];

   delete[] a;
   delete[] b;
   delete[] part_c;

   return dot_product;
}

const int coutWidth = 14;
const int numOfTries = 5;

int main()
{
   std::cout << std::scientific;
   std::cout << std::setw(coutWidth) << "true dot" << std::setw(coutWidth) << "GPU computed" << std::setw(coutWidth) << "diff" << std::setw(coutWidth) << "time" << std::endl;

   //double true_dot = 2 * sumOfSquares(N - 1);
   double true_dot = sumOfN(N);
   double mean_time = 0;

   for(size_t i = 0; i < numOfTries; i++)
   {
      float time = 0;
      double gpu_dot = testDot(time);
      std::cout << std::setw(coutWidth) << true_dot << std::setw(coutWidth) << gpu_dot << std::setw(coutWidth) << abs(gpu_dot - true_dot) << std::setw(coutWidth) << time << std::endl;
      mean_time += time;
   }

   std::cout << "Mean time is: " << mean_time / numOfTries << std::endl;

   return 0;
}