#include "laplace_solver.h"

const double PI = 3.14159265359;

double laplace_solver::u(Vector3 v)
{
   return v.x * v.x - v.y * v.y;
}

Vector3 laplace_solver::gradU(Vector3 v)
{
   return { 2 * v.x, -2 * v.y, 0 };
}

double laplace_solver::laplaceIntegral1(Vector3 v,
                                       Vector3 point,
                                       Vector3 normal)
{
   Vector3 grad = gradU(v);
   double dudnx = grad.x * normal.x;
   double dudny = grad.y * normal.y;
   double dudnz = grad.z * normal.z;

   return 1 / (v - point).GetLength() * (dudnx + dudny + dudnz);
}

double laplace_solver::laplaceIntegral2(Vector3 v,
                                        Vector3 point,
                                        Vector3 normal)
{
   double l = (v - point).GetLength();

   double rx = normal.x / ((v.x - point.x) * l);
   double ry = normal.y / ((v.y - point.y) * l);
   double rz = normal.z / ((v.z - point.z) * l);

   return -1 * (rx + ry + rz) * u(v);
}

void laplace_solver::calcIntegralOverMesh(const Mesh& mesh,
                                         const QuadPoints& qp,
                                         const vector<Vector3>& points,
                                         vector<double>& result)
{
   result = vector<double>(points.size());

   for (size_t p = 0; p < points.size(); p++)
   {
      double integral_1 = 0;
      double integral_2 = 0;

      for (size_t t = 0; t < mesh.TriangleCount(); t++)
      {
         double tringle_sum_1 = 0;
         double tringle_sum_2 = 0;

         Triangle tr = mesh.GetTriangle(t);

         for (size_t o = 0; o < qp.order; o++)
         {
            Vector3 v = tr.PointFromST(qp.x[o], qp.y[o]);
            tringle_sum_1 += qp.w[o] * laplaceIntegral1(v, points[p], tr.Normal());
            tringle_sum_2 += qp.w[o] * laplaceIntegral2(v, points[p], tr.Normal());
         }

         integral_1 += tringle_sum_1 * tr.Area();
         integral_2 += tringle_sum_2 * tr.Area();
      }

      result[p] = 1.0 / (4.0 * PI) * (integral_1 - integral_2);
   }
}
