#include "laplace_solver.h"

const double PI = 3.14159265359;

double laplace_solver::u(Vector3 v)
{
   return 2 * v.x * v.x - v.y * v.y - v.z * v.z;
   //return 1;
   //return v.x;
}

Vector3 laplace_solver::gradU(Vector3 v)
{
   return { 4 * v.x, -2 * v.y, -2 * v.z };
   //return { 1.0, 0.0, 0.0 };
}

double laplace_solver::laplaceIntegral1(Vector3 v,
                                        Vector3 point,
                                        Vector3 normal)
{
   Vector3 grad = gradU(v);

   double dudnx = grad.x * normal.x;
   double dudny = grad.y * normal.y;
   double dudnz = grad.z * normal.z;

   return (dudnx + dudny + dudnz) / (point - v).Length();
}

double laplace_solver::laplaceIntegral2(Vector3 v,
                                        Vector3 point,
                                        Vector3 normal)
{
   double l = (point - v).Length();

   double rx = normal.x * (point.x - v.x);
   double ry = normal.y * (point.y - v.y);
   double rz = normal.z * (point.z - v.z);

   return (rx + ry + rz) * u(v) / pow(l, 3.0);
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
         Vector3 normal = tr.Normal();

         for (size_t o = 0; o < qp.order; o++)
         {
            Vector3 v = tr.PointFromST(qp.x[o], qp.y[o]);
            tringle_sum_1 += qp.w[o] * laplaceIntegral1(v, points[p], normal);
            tringle_sum_2 += qp.w[o] * laplaceIntegral2(v, points[p], normal);
         }

         integral_1 += tringle_sum_1 * tr.Area();
         integral_2 += tringle_sum_2 * tr.Area();
      }

      result[p] = (integral_1 - integral_2) / (4.0 * PI);
   }
}
