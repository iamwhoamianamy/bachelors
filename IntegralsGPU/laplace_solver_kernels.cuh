#pragma once
#include "real.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "vector3.cuh"
#include "laplace_solver_structs.h"

using namespace triangle_quadratures;

namespace laplace_solver_kernels
{
   __global__ void solverKernelArraysReduction(
      const real* quadratures_X,
      const real* quadratures_Y,
      const real* quadratures_Z,
      const real* normals_X,
      const real* normals_Y,
      const real* normals_Z,
      const real* points_X,
      const real* points_Y,
      const real* points_Z,
      const real* weights,
      const real* areas,
      const int trianglesCount,
      const int pointsCount,
      const int quadraturesOrder,
      real* results);

   __global__ void solverKernelVector3sReduction(
      const Vector3* quadPoints,
      const Vector3* normals,
      const Vector3* points,
      const real* weights,
      const real* areas,
      const int trianglesCount,
      const int pointsCount,
      const int quadraturesOrder,
      real* results);

   __global__ void solverKernelVector3sBlocks(
      const Vector3* quadPoints,
      const Vector3* normals,
      const Vector3* points,
      const real* weights,
      const real* areas,
      const int trianglesCount,
      const int pointsCount,
      const int quadraturesOrder,
      real* results);

   __global__ void solverKernelStructsReduction(
      const QuadPoint* quadPoints,
      const Vector3* points,
      const int quadraturesCount,
      real* results);

   __global__ void solverKernelStructsBlocks(
      const QuadPoint* quadPoints,
      const Vector3* points,
      const int quadraturesCount,
      real* results);

   __global__ void solverKernelStructsGrid(
      const QuadPoint* quadPoints,
      const Vector3* points,
      const int matrixWidth,
      real* resultsMatrix);

   //__global__ void AddMatrices(const real* a, const real* b, real* c);
   //__global__ void AddMatricesShared(const double* a, const double* b, double* c);

}