#include "laplace_solver.h"

const float PI = 3.14159265359;

float laplace_solver::u(Vector3 v)
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

float laplace_solver::laplaceIntegral1(Vector3 v,
                                        Vector3 point,
                                        Vector3 normal)
{
   Vector3 grad = gradU(v);

   float dudnx = grad.x * normal.x;
   float dudny = grad.y * normal.y;
   float dudnz = grad.z * normal.z;

   return (dudnx + dudny + dudnz) / (point - v).Length();
}

float laplace_solver::laplaceIntegral2(Vector3 v,
                                        Vector3 point,
                                        Vector3 normal)
{
   float l = (point - v).Length();

   float rx = normal.x * (point.x - v.x);
   float ry = normal.y * (point.y - v.y);
   float rz = normal.z * (point.z - v.z);

   return (rx + ry + rz) * u(v) / pow(l, 3.0);
}

void laplace_solver::calcIntegralOverMesh(const Mesh& mesh,
                                          const BasisQuadratures& qp,
                                          const vector<Vector3>& points,
                                          vector<float>& results)
{
   results = vector<float>(points.size());

   for (size_t p = 0; p < points.size(); p++)
   {
      float integral1 = 0;
      float integral_2 = 0;

      for (size_t t = 0; t < mesh.TrianglesCount(); t++)
      {
         float tringle_sum_1 = 0;
         float tringle_sum_2 = 0;

         Triangle tr = mesh.GetTriangle(t);
         Vector3 normal = tr.Normal();

         for (size_t o = 0; o < qp.order; o++)
         {
            Vector3 v = tr.PointFromST(qp.x[o], qp.y[o]);
            tringle_sum_1 += qp.w[o] * laplaceIntegral1(v, points[p], normal);
            tringle_sum_2 += qp.w[o] * laplaceIntegral2(v, points[p], normal);
         }

         integral1 += tringle_sum_1 * tr.Area();
         integral_2 += tringle_sum_2 * tr.Area();
      }

      results[p] = (integral1 - integral_2) / (4.0 * PI);
   }
}
