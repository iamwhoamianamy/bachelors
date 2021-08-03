#include "bitmap_image.hpp"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "fractal.cuh"

__device__ void SetPixel(unsigned char* ptr, const int offset, const unsigned char r, const unsigned char g, const unsigned char b)
{
   ptr[offset * 3 + 0] = r;
   ptr[offset * 3 + 1] = g;
   ptr[offset * 3 + 2] = b;
}

__device__ int Julia(const int x, const int y)
{
   const float SCALE = 1.5;

   float jx = SCALE * (float)(DIM / 2 - x) / (DIM / 2);
   float jy = SCALE * (float)(DIM / 2 - y) / (DIM / 2);

   Complex c(-0.8, 0.156);
   Complex a(jx, jy);

   size_t i;
   for(i = 0; i < 200; i++)
   {
      a = a * a + c;
      if(a.MagSquared() > 1000)
         return 0;
   }

   return 1;
}

__global__ void kernel(unsigned char* ptr)
{
   int x = blockIdx.x;
   int y = blockIdx.y;
   int offset = x + y * gridDim.x;

   int julia_value = Julia(x, y);
   SetPixel(ptr, offset, 255 * julia_value, 0, 0);
}

__device__ void FormFractal()
{
   bitmap_image image(DIM, DIM);
   image.clear();

   unsigned char* dev_bitmap;
   unsigned char* bitmap_pixels = new unsigned char[IMAGE_SIZE];

   TryCudaMalloc((void**)&dev_bitmap, IMAGE_SIZE);

   dim3 grid(DIM, DIM);

   kernel<<<grid, 1>>>(dev_bitmap);

   HANDLE_ERROR(cudaMemcpy(bitmap_pixels, dev_bitmap, IMAGE_SIZE, cudaMemcpyDeviceToHost));
   cudaFree(dev_bitmap);

   for(size_t y = 0; y < DIM; y++)
   {
      for(size_t x = 0; x < DIM; x++)
      {
         int offset = x + y * DIM;
         image.set_pixel(x, y, bitmap_pixels[offset * 3 + 0], 0, 0);
      }
   }

   delete[] bitmap_pixels;

   image.save_image("image.bmp");
}