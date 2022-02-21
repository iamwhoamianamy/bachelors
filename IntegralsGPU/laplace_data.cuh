#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "vector3.cuh"
#include "real.h"

using namespace triangle_quadratures;

namespace laplace_data
{
   __device__ const real PI = 3.14159265359;
   __device__ const real fourPI = 4 * 3.14159265359;

   __device__ __host__ inline real u(
      const real x, const real y, const real z);

   __device__ __host__ inline real gradUX(
      const real x, const real y, const real z);

   __device__ __host__ inline real gradUY(
      const real x, const real y, const real z);

   __device__ __host__ inline real gradUZ(
      const real x, const real y, const real z);

   __device__ inline real calcInverseDistanceGPU(
      const real x1, const real y1, const real z1,
      const real x2, const real y2, const real z2);

   __host__ inline real calcDistanceCPU(
      const real x1, const real y1, const real z1,
      const real x2, const real y2, const real z2);

   __device__ real calcLaplaceIntegralGPU(
      const real qx, const real qy, const real qz,
      const real px, const real py, const real pz,
      const real nx, const real ny, const real nz);

   __host__ real calcLaplaceIntegralCPU(
      const real qx, const real qy, const real qz,
      const real px, const real py, const real pz,
      const real nx, const real ny, const real nz);

   const int MAX_THREADS = 1024;

   const int THREADS_PER_BLOCK = 512;

   const int QUADS_PER_BLOCK = 32;
   const int QUAD_ORDER = 64;
   //const int TR_PER_BLOCK = QUADS_PER_BLOCK / QUAD_ORDER;
   const int TR_PER_BLOCK = 1;

   const int POINTS_PER_BLOCK = MAX_THREADS / THREADS_PER_BLOCK;

   const int BLOCK_SIZE = 128;
};