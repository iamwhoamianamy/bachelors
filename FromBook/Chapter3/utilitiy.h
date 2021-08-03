#include "cudahelper.h"

void SetBestFitDevice()
{
   cudaDeviceProp info;
   int dev_index;

   HANDLE_ERROR(cudaGetDevice(&dev_index));
   printf("Current device index: %d\n", dev_index);

   info = cudaDeviceProp();
   info.major = 1;
   info.minor = 3;
   HANDLE_ERROR(cudaChooseDevice(&dev_index, &info));
   printf("The closest device index to 1:3 level is %d\n", dev_index);
   HANDLE_ERROR(cudaSetDevice(dev_index));
}

__global__ void add(int a, int b, int* c)
{
   *c = a + b;
}

void AddTwoNumbers()
{
   int c;
   int* dev_c;
   HANDLE_ERROR(cudaMalloc((void**)&dev_c, sizeof(int)));

   add<<<1,1>>>(2, 7, dev_c);

   HANDLE_ERROR(cudaMemcpy(&c, dev_c, sizeof(int), cudaMemcpyDeviceToHost));

   printf("2 + 7 = %d\n", c);
   cudaFree(dev_c);
}