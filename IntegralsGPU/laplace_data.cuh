#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "vector3.cuh"

using namespace triangle_quadratures;

#ifndef LAPLACE_DATA
#define LAPLACE_DATA

namespace laplace_data
{
   __device__ const double PI = 3.14159265359;

   __device__ __host__ double u(double x, double y, double z);
   __device__ __host__ double gradUX(double x, double y, double z);
   __device__ __host__ double gradUY(double x, double y, double z);
   __device__ __host__ double gradUZ(double x, double y, double z);

   __device__ __host__ double lengthBetween(double x1, double y1, double z1,
                                            double x2, double y2, double z2);

   __device__ __host__ double laplaceIntegral1(double qx, double qy, double qz,
                                               double px, double py, double pz,
                                               double nx, double ny, double nz);

   __device__ __host__ double laplaceIntegral2(double qx, double qy, double qz,
                                               double px, double py, double pz,
                                               double nx, double ny, double nz);

   __device__ __host__ double laplaceIntegral1(const Vector3& q,
                                               const Vector3& p,
                                               const Vector3& n);

   __device__ __host__ double laplaceIntegral2(const Vector3& q,
                                               const Vector3& p,
                                               const Vector3& n);
};

#endif // !LAPLACE_DATA