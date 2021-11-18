#include "laplace_data.cuh"
#include <cmath>


__device__ __host__ float laplace_data::u(float x, float y, float z)
{
   return 2 * x * x - y * y - z * z;
}

__device__ __host__ float laplace_data::gradUX(float x, float y, float z)
{
   return 4 * x;
}

__device__ __host__ float laplace_data::gradUY(float x, float y, float z)
{
   return -2 * y;
}

__device__ __host__ float laplace_data::gradUZ(float x, float y, float z)
{
   return -2 * z;
}

__device__ __host__ float laplace_data::lengthBetween(float x1, float y1, float z1,
                                                      float x2, float y2, float z2)
{
   return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
}

__device__ __host__ float laplace_data::laplaceIntegral1(float qx, float qy, float qz,
                                                         float px, float py, float pz,
                                                         float nx, float ny, float nz)
{
   float gradx = gradUX(qx, qy, qz);
   float grady = gradUY(qx, qy, qz);
   float gradz = gradUZ(qx, qy, qz);

   float dudnx = gradx * nx;
   float dudny = grady * ny;
   float dudnz = gradz * nz;

   return (dudnx + dudny + dudnz) / lengthBetween(qx, qy, qz, px, py, pz);
}

__device__ __host__ float laplace_data::laplaceIntegral2(float qx, float qy, float qz,
                                                         float px, float py, float pz,
                                                         float nx, float ny, float nz)
{
   float l = lengthBetween(qx, qy, qz, px, py, pz);

   float rx = nx * (px - qx);
   float ry = ny * (py - qy);
   float rz = nz * (pz - qz);

   return (rx + ry + rz) * u(qx, qy, qz) / (l * l * l);
}

/*
__device__ __host__ float laplace_data::u(Vector3 vec)
{
   return u(vec.x, vec.y, vec.z);
}

__device__ __host__ Vector3 laplace_data::gradU(Vector3 vec)
{
   return Vector3(gradUX(vec.x, vec.y, vec.z),
                  gradUY(vec.x, vec.y, vec.z),
                  gradUZ(vec.x, vec.y, vec.z));
}*/

__device__ __host__ float laplace_data::laplaceIntegral1(const Vector3& q, const Vector3& p, const Vector3& n)
{
   return laplaceIntegral1(q.x, q.y, q.z, p.x, p.y, p.z, n.x, n.y, n.z);
}

__device__ __host__ float laplace_data::laplaceIntegral2(const Vector3& q, const Vector3& p, const Vector3& n)
{
   return laplaceIntegral2(q.x, q.y, q.z, p.x, p.y, p.z, n.x, n.y, n.z);
}