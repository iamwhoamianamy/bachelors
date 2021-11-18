#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "vector3.cuh"

using namespace triangle_quadratures;

namespace laplace_data
{
   __device__ const float PI = 3.14159265359;

   __device__ __host__ inline float u(float x, float y, float z);
   __device__ __host__ inline float gradUX(float x, float y, float z);
   __device__ __host__ inline float gradUY(float x, float y, float z);
   __device__ __host__ inline float gradUZ(float x, float y, float z);

   __device__ __host__ float lengthBetween(float x1, float y1, float z1,
                                            float x2, float y2, float z2);

   __device__ __host__ float laplaceIntegral1(float qx, float qy, float qz,
                                               float px, float py, float pz,
                                               float nx, float ny, float nz);

   __device__ __host__ float laplaceIntegral2(float qx, float qy, float qz,
                                               float px, float py, float pz,
                                               float nx, float ny, float nz);

   __device__ __host__ float laplaceIntegral1(const Vector3& q,
                                               const Vector3& p,
                                               const Vector3& n);

   __device__ __host__ float laplaceIntegral2(const Vector3& q,
                                               const Vector3& p,
                                               const Vector3& n);

   const int MAX_THREADS = 1024;

   const int THREADS_PER_BLOCK = 512;

   const int QUADS_PER_BLOCK = 32;
   const int QUAD_ORDER = 64;
   //const int TR_PER_BLOCK = QUADS_PER_BLOCK / QUAD_ORDER;
   const int TR_PER_BLOCK = 1;

   const int POINTS_PER_BLOCK = MAX_THREADS / THREADS_PER_BLOCK;
};