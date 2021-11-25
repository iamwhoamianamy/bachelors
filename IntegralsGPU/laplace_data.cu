#include "laplace_data.cuh"
#include <cmath>


__device__ __host__ real laplace_data::u(const real x, const real y, const real z)
{
   return 2 * x * x - y * y - z * z;
}

__device__ __host__ real laplace_data::gradUX(const real x, const real y, const real z)
{
   return 4 * x;
}

__device__ __host__ real laplace_data::gradUY(const real x, const real y, const real z)
{
   return -2 * y;
}

__device__ __host__ real laplace_data::gradUZ(const real x, const real y, const real z)
{
   return -2 * z;
}

__device__ __host__ real laplace_data::lengthBetween(const real x1, const real y1, const real z1,
                                                     const real x2, const real y2, const real z2)
{
   return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
}

__device__ __host__ real laplace_data::laplaceIntegral1(const real qx, const real qy, const real qz,
                                                        const real px, const real py, const real pz,
                                                        const real nx, const real ny, const real nz)
{
   real dudnx = gradUX(qx, qy, qz) * nx;
   real dudny = gradUY(qx, qy, qz) * ny;
   real dudnz = gradUZ(qx, qy, qz) * nz;

   return (dudnx + dudny + dudnz) / lengthBetween(qx, qy, qz, px, py, pz);
}

__device__ __host__ real laplace_data::laplaceIntegral2(const real qx, const real qy, const real qz,
                                                        const real px, const real py, const real pz,
                                                        const real nx, const real ny, const real nz)
{
   real l = lengthBetween(qx, qy, qz, px, py, pz);

   real rx = nx * (px - qx);
   real ry = ny * (py - qy);
   real rz = nz * (pz - qz);

   return (rx + ry + rz) * u(qx, qy, qz) / (l * l * l);
}

/*
__device__ __host__ real laplace_data::u(Vector3 vec)
{
   return u(vec.x, vec.y, vec.z);
}

__device__ __host__ Vector3 laplace_data::gradU(Vector3 vec)
{
   return Vector3(gradUX(vec.x, vec.y, vec.z),
                  gradUY(vec.x, vec.y, vec.z),
                  gradUZ(vec.x, vec.y, vec.z));
}*/

__device__ __host__ real laplace_data::laplaceIntegral1(const Vector3& q, const Vector3& p, const Vector3& n)
{
   return laplaceIntegral1(q.x, q.y, q.z, p.x, p.y, p.z, n.x, n.y, n.z);
}

__device__ __host__ real laplace_data::laplaceIntegral2(const Vector3& q, const Vector3& p, const Vector3& n)
{
   return laplaceIntegral2(q.x, q.y, q.z, p.x, p.y, p.z, n.x, n.y, n.z);
}