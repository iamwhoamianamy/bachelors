#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cudahelper.cuh"

#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <chrono>

//const int SIZE = 1024;

float calcHistoOnCPU(const std::string& text, unsigned int* histo)
{
   auto start = std::chrono::steady_clock::now();

   for(auto c : text)
      histo[c]++;

   auto stop = std::chrono::steady_clock::now();

   return (stop - start).count() * 1e-12;
}

__global__ void histoKernel(unsigned char* text, int size, unsigned int* histo)
{
   __shared__ unsigned int temp[256];

   temp[threadIdx.x] = 0;
   __syncthreads();

   int i = threadIdx.x + blockIdx.x * blockDim.x;
   int offset = blockDim.x * gridDim.x;

   while(i < size)
   {
      atomicAdd(&temp[text[i]], 1);
      i += offset;
   }

   __syncthreads();

   atomicAdd(&histo[threadIdx.x], temp[threadIdx.x]);
}

float calcHistoOnGPU(const std::string& text, unsigned int* histo)
{
   const int text_size = text.size();

   unsigned char* char_text = (unsigned char*)text.begin()._Ptr;
   unsigned char* dev_char_text;
   
   tryCudaMalloc((void**)&dev_char_text, text_size * sizeof(char));
   handleError(cudaMemcpy(dev_char_text, char_text, text_size * sizeof(char), cudaMemcpyHostToDevice));

   unsigned int* dev_histo;

   tryCudaMalloc((void**)&dev_histo, 256 * sizeof(int));
   handleError(cudaMemset(dev_histo, 0, 256 * sizeof(int)));

   cudaEvent_t start, stop;
   handleError(cudaEventCreate(&start));
   handleError(cudaEventCreate(&stop));
   handleError(cudaEventRecord(start, 0));

   histoKernel<<<64, 256>>>(dev_char_text, text_size, dev_histo);
   tryKernelLaunch();
   tryKernelSynchronize();

   handleError(cudaEventRecord(stop, 0));
   handleError(cudaEventSynchronize(stop));

   float elapsed_time = 0;
   handleError(cudaEventElapsedTime(&elapsed_time, start, stop));

   handleError(cudaEventDestroy(start));
   handleError(cudaEventDestroy(stop));

   handleError(cudaMemcpy(histo, dev_histo, 256 * sizeof(int), cudaMemcpyDeviceToHost));

   cudaFree(dev_histo);
   cudaFree(dev_char_text);

   return elapsed_time * 1e-6;
}

std::string fileToString(const std::string& fileName)
{
   std::ifstream fin;
   fin.open(fileName, std::ifstream::in);

   std::string content((std::istreambuf_iterator<char>(fin)),
                        std::istreambuf_iterator<char>());

   fin.close();

   return content;
}

void printHistograms(const unsigned int* histoCPU, const unsigned int* histoGPU)
{
   for(size_t i = 0; i < 256; i++)
   {
      std::cout << char(i) << " - " << histoCPU[i] << " " << histoGPU[i] << std::endl;
   }
}

int main()
{
   auto str = fileToString("text.txt");

   unsigned int* histo_CPU = new unsigned int[256] { 0 };
   float elapsed_CPU = calcHistoOnCPU(str, histo_CPU);

   unsigned int* histo_GPU = new unsigned int[256] { 0 };
   float elapsed_GPU = calcHistoOnGPU(str, histo_GPU);

   printHistograms(histo_CPU, histo_GPU);

   std::cout << "CPU time: " << std::scientific << elapsed_CPU << std::fixed << std::endl;
   std::cout << "GPU time: " << std::scientific << elapsed_GPU << std::fixed << std::endl;

   std::cout << "Acceleration with GPU is " << elapsed_CPU / elapsed_GPU;

   return 0;
}