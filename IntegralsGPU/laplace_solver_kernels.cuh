#pragma once
#include "real.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "vector3.cuh"
#include "laplace_solver_structs.h"

using namespace triangle_quadratures;

namespace laplace_solver_kernels
{
   __global__ void solverKernelVector3s(
      const Vector3* quadPoints,
      const Vector3* normals,
      const Vector3* points,
      const real* weights,
      const real* areas,
      const int trianglesCount,
      const int pointsCount,
      const int quadraturesOrder,
      real* results);

   __global__ void solverKernelStructs(
      const QuadPoint* quadPoints,
      const Vector3* points,
      const int quadraturesCount,
      real* results);
}