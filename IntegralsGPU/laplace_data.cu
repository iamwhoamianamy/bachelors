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
