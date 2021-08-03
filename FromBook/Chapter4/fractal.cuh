#include "bitmap_image.hpp"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cudahelper.h"

const int DIM = 500;
const int IMAGE_SIZE = DIM * DIM * 3;

struct Complex
{
   float r;
   float i;

   __device__ Complex(float a, float b) : r(a), i(b) {}

   __device__ float MagSquared()
   {
      return r * r + i * i;
   }

   __device__ Complex operator * (const Complex& rhs)
   {
      return Complex(r * rhs.r - i * rhs.i, i * rhs.r + r * rhs.i);
   }

   __device__ Complex operator + (const Complex& rhs)
   {
      return Complex(r + rhs.r, i + rhs.i);
   }
};

__device__ void SetPixel(unsigned char* ptr, const int offset, const unsigned char r, const unsigned char g, const unsigned char b);

__device__ int Julia(const int x, const int y);

__global__ void kernel(unsigned char* ptr);

__device__ void FormFractal();