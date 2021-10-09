#include "laplace_solver.h"

double LaplaceSolver::U(Vector3 v)
{
   return v.x * v.x - v.y * v.y;
}

Vector3 LaplaceSolver::GradU(Vector3 v)
{
   return { 2 * v.x, -2 * v.y, 0 };
}

double LaplaceSolver::LaplaceIntegral1(Vector3 v,
                                       Vector3 point,
                                       Vector3 normal)
{
   Vector3 grad = GradU(v);
   double dudnx = grad.x * normal.x;
   double dudny = grad.y * normal.y;
   double dudnz = grad.z * normal.z;

   return 1 / (v - point).GetLength() * (dudnx + dudny + dudnz);
}

double LaplaceSolver::LaplaceIntegral2(Vector3 v,
                                       Vector3 point,
                                       Vector3 normal)
{
   double l = (v - point).GetLength();

   double rx = normal.x / ((v.x - point.x) * l);
   double ry = normal.y / ((v.y - point.y) * l);
   double rz = normal.z / ((v.z - point.z) * l);

   return (rx + ry + rz) * U(v);
}
