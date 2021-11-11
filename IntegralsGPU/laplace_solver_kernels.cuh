#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "vector3.cuh"

using namespace triangle_quadratures;

namespace laplace_solver_kernels
{
   __global__ void SolverKernelArrays(
      const double* quadratures_X,
      const double* quadratures_Y,
      const double* quadratures_Z,
      const double* normals_X,
      const double* normals_Y,
      const double* normals_Z,
      const double* points_X,
      const double* points_Y,
      const double* points_Z,
      const double* weights,
      const double* areas,
      const int trianglesCount,
      const int pointsCount,
      const int quadraturesOrder,
      double* result);

   __global__ void SolverKernelVector3s(
      const Vector3* quadratures,
      const Vector3* normals,
      const Vector3* points,
      const double* weights,
      const double* areas,
      const int trianglesCount,
      const int pointsCount,
      const int quadraturesOrder,
      double* result);
}