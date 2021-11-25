#include "laplace_data.cuh"
#include <cmath>


__device__ __host__ inline real laplace_data::u(const real x, const real y, const real z)
{
   return 2 * x * x - y * y - z * z;
}

__device__ __host__ inline real laplace_data::gradUX(const real x, const real y, const real z)
{
   return 4 * x;
}

__device__ __host__ inline real laplace_data::gradUY(const real x, const real y, const real z)
{
   return -2 * y;
}

__device__ __host__ inline real laplace_data::gradUZ(const real x, const real y, const real z)
{
   return -2 * z;
}

__device__ float rsqrtf(float  x);
__device__ float rsqrtf(float  x);
__device__ double rsqrt(double  x);
__device__ float norm3df(float  a, float  b, float  c);
__device__ float rnorm3df(float  a, float  b, float  c);

__device__ real laplace_data::inverseDistanceGPU(const real x1, const real y1, const real z1,
                                               const real x2, const real y2, const real z2)
{
#ifdef REAL_IS_FLOAT
   return rsqrtf((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
   //return sqrtf((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
   //return norm3df(x1 - x2, y1 - y2, z1 - z2);
   //return rnorm3df(x1 - x2, y1 - y2, z1 - z2);
#endif // REAL_IS_FLOAT

#ifdef REAL_IS_DOUBLE
   return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
#endif // REAL_IS_DOUBLE
}

__device__ real laplace_data::laplaceIntegral1GPU(const real qx, const real qy, const real qz,
                                                  const real px, const real py, const real pz,
                                                  const real nx, const real ny, const real nz)
{
   real dudnx = gradUX(qx, qy, qz) * nx;
   real dudny = gradUY(qx, qy, qz) * ny;
   real dudnz = gradUZ(qx, qy, qz) * nz;

   real l = inverseDistanceGPU(qx, qy, qz, px, py, pz);

   return (dudnx + dudny + dudnz) * l;
}

__device__ real laplace_data::laplaceIntegral2GPU(const real qx, const real qy, const real qz,
                                                  const real px, const real py, const real pz,
                                                  const real nx, const real ny, const real nz)
{
   real l = inverseDistanceGPU(qx, qy, qz, px, py, pz);

   real rx = nx * (px - qx);
   real ry = ny * (py - qy);
   real rz = nz * (pz - qz);

   return (rx + ry + rz) * u(qx, qy, qz) * (l * l * l);
}

__device__ real laplace_data::laplaceIntegral1GPU(const Vector3& q, const Vector3& p, const Vector3& n)
{
   return laplaceIntegral1GPU(q.x, q.y, q.z, p.x, p.y, p.z, n.x, n.y, n.z);
}

__device__ real laplace_data::laplaceIntegral2GPU(const Vector3& q, const Vector3& p, const Vector3& n)
{
   return laplaceIntegral2GPU(q.x, q.y, q.z, p.x, p.y, p.z, n.x, n.y, n.z);
}

__host__ real laplace_data::lengthBetweenCPU(const real x1, const real y1, const real z1,
                                             const real x2, const real y2, const real z2)
{
#ifdef REAL_IS_FLOAT
   return sqrtf((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
#endif // REAL_IS_FLOAT

#ifdef REAL_IS_DOUBLE
   return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
#endif // REAL_IS_DOUBLE
}

__host__ real laplace_data::laplaceIntegral1CPU(const real qx, const real qy, const real qz,
                                                const real px, const real py, const real pz,
                                                const real nx, const real ny, const real nz)
{
   real dudnx = gradUX(qx, qy, qz) * nx;
   real dudny = gradUY(qx, qy, qz) * ny;
   real dudnz = gradUZ(qx, qy, qz) * nz;

   return (dudnx + dudny + dudnz) / lengthBetweenCPU(qx, qy, qz, px, py, pz);
}

__host__ real laplace_data::laplaceIntegral2CPU(const real qx, const real qy, const real qz,
                                                const real px, const real py, const real pz,
                                                const real nx, const real ny, const real nz)
{
   real l = lengthBetweenCPU(qx, qy, qz, px, py, pz);

   real rx = nx * (px - qx);
   real ry = ny * (py - qy);
   real rz = nz * (pz - qz);

   return (rx + ry + rz) * u(qx, qy, qz) / (l * l * l);
}

__host__ real laplace_data::laplaceIntegral1CPU(const Vector3& q, const Vector3& p, const Vector3& n)
{
   return laplaceIntegral1CPU(q.x, q.y, q.z, p.x, p.y, p.z, n.x, n.y, n.z);
}

__host__ real laplace_data::laplaceIntegral2CPU(const Vector3& q, const Vector3& p, const Vector3& n)
{
   return laplaceIntegral2CPU(q.x, q.y, q.z, p.x, p.y, p.z, n.x, n.y, n.z);
}