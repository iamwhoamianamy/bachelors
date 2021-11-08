#include "laplace_data.cuh"
#include <cmath>


__device__ double laplace_data::u(double x, double y, double z)
{
   return 2 * x * x - y * y - z * z;
}

__device__ double laplace_data::gradUX(double x, double y, double z)
{
   return 4 * x;
}

__device__ double laplace_data::gradUY(double x, double y, double z)
{
   return -2 * y;
}

__device__ double laplace_data::gradUZ(double x, double y, double z)
{
   return -2 * z;
}

__device__ double laplace_data::lengthBetween(double x1, double y1, double z1,
                                              double x2, double y2, double z2)
{
   return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
}

__device__ double laplace_data::laplaceIntegral1(double qx, double qy, double qz,
                                                 double px, double py, double pz,
                                                 double nx, double ny, double nz)
{
   double gradx = gradUX(qx, qy, qz);
   double grady = gradUY(qx, qy, qz);
   double gradz = gradUZ(qx, qy, qz);

   double dudnx = gradx * nx;
   double dudny = grady * ny;
   double dudnz = gradz * nz;

   return (dudnx + dudny + dudnz) / lengthBetween(qx, qy, qz, px, py, pz);
}

__device__ double laplace_data::laplaceIntegral2(double qx, double qy, double qz,
                                                 double px, double py, double pz,
                                                 double nx, double ny, double nz)
{
   double l = lengthBetween(qx, qy, qz, px, py, pz);

   double rx = nx * (px - qx);
   double ry = ny * (py - qy);
   double rz = nz * (pz - qz);

   return (rx + ry + rz) * u(qx, qy, qz) / (l * l * l);
}

/*
__device__ __host__ double laplace_data::u(Vector3 vec)
{
   return u(vec.x, vec.y, vec.z);
}

__device__ __host__ Vector3 laplace_data::gradU(Vector3 vec)
{
   return Vector3(gradUX(vec.x, vec.y, vec.z),
                  gradUY(vec.x, vec.y, vec.z),
                  gradUZ(vec.x, vec.y, vec.z));
}*/

__device__ __host__ double laplace_data::laplaceIntegral1(const Vector3& q, const Vector3& p, const Vector3& n)
{
   return laplaceIntegral1(q.x, q.y, q.z, p.x, p.y, p.z, n.x, n.y, n.z);
}

__device__ __host__ double laplace_data::laplaceIntegral2(const Vector3& q, const Vector3& p, const Vector3& n)
{
   return laplaceIntegral2(q.x, q.y, q.z, p.x, p.y, p.z, n.x, n.y, n.z);
}

