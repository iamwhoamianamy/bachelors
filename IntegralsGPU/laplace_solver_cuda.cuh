#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <vector>
#include "vector3.h"
#include "mesh.h"
#include "quad_points.h"

using namespace triangle_quadratures;
using namespace std;

namespace laplace_solver_cuda
{
   __device__ double u(double x, double y, double z);
   __device__ double gradUX(double x, double y, double z);
   __device__ double gradUY(double x, double y, double z);
   __device__ double gradUZ(double x, double y, double z);

   __device__ double lengthBetween(double x1, double y1, double z1,
                                   double x2, double y2, double z2);

   __device__ double laplaceIntegral1(double qx, double qy, double qz,
                                      double px, double py, double pz,
                                      double nx, double ny, double nz);

   __device__ double laplaceIntegral2(double qx, double qy, double qz,
                                      double px, double py, double pz,
                                      double nx, double ny, double nz);

   void calcIntegralOverMesh(const Mesh& mesh,
                             const QuadPoints& qp,
                             const vector<Vector3>& points,
                             vector<double>& result);

   __global__ void calcIntegralOverMeshArrays(const double* quadraturesX,
                                   const double* quadraturesY,
                                   const double* quadraturesZ,
                                   const double* normalsX,
                                   const double* normalsY,
                                   const double* normalsZ,
                                   const double* pointsX,
                                   const double* pointsY, 
                                   const double* pointsZ,
                                   const double* weights,
                                   const double* areas,
                                   const int trianglesCount,
                                   const int pointsCount,
                                   const int quadratureOrder,
                                   double* result);

   typedef unsigned int uint;
};