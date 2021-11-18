#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "vector3.cuh"
#include "laplace_solver_structs.h"

using namespace triangle_quadratures;

namespace laplace_solver_kernels
{
   __global__ void SolverKernelArraysReduction(
      const float* quadratures_X,
      const float* quadratures_Y,
      const float* quadratures_Z,
      const float* normals_X,
      const float* normals_Y,
      const float* normals_Z,
      const float* points_X,
      const float* points_Y,
      const float* points_Z,
      const float* weights,
      const float* areas,
      const int trianglesCount,
      const int pointsCount,
      const int quadraturesOrder,
      float* results);

   __global__ void SolverKernelVector3sReduction(
      const Vector3* quadPoints,
      const Vector3* normals,
      const Vector3* points,
      const float* weights,
      const float* areas,
      const int trianglesCount,
      const int pointsCount,
      const int quadraturesOrder,
      float* results);

   __global__ void SolverKernelVector3sBlocks(
      const Vector3* quadPoints,
      const Vector3* normals,
      const Vector3* points,
      const float* weights,
      const float* areas,
      const int trianglesCount,
      const int pointsCount,
      const int quadraturesOrder,
      float* results);

   __global__ void SolverKernelStructsReduction(
      const QuadPoint* quadPoints,
      const Vector3* points,
      const int quadraturesCount,
      float* results);

   __global__ void SolverKernelStructsBlocks(
      const QuadPoint* quadPoints,
      const Vector3* points,
      const int trianglesCount,
      const int pointsCount,
      const int quadraturesOrder,
      float* results);
}