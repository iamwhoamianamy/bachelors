#include "triangle_quadratures.h"
#include <vector>

namespace triangle_quadratures
{
   double calcIntegralOverTriangle(double (*f)(double, double, double, double),
                                   const Triangle& tr,
                                   const QuadPoints& qp)
   {
      double result = 0;

      for (size_t i = 0; i < qp.order; i++)
      {
         Point p = tr.PointFromST(qp.x[i], qp.y[i]);
         result += qp.w[i] * f(p.x, p.y, p.z, 0);
      }

      return result * tr.Area();
   }

   void calcIntegralOverMesh(double(*f)(double, double, double, double),
                             const Mesh& mesh, 
                             const QuadPoints& qp,
                             const vector<Point>& points,
                             std::vector<double>& result)
   {
      vector<double> points_double(points.size() * 3);

      for(size_t i = 0; i < points.size(); i++)
      {
         points_double[i * 3 + 0] = points[i].x;
         points_double[i * 3 + 1] = points[i].y;
         points_double[i * 3 + 2] = points[i].z;
      }

      vector<double> coords(mesh.TriangleCount() * qp.order * 3);
      vector<double> areas(mesh.TriangleCount());

      for(size_t i = 0; i < mesh.TriangleCount(); i++)
      {
         Triangle tr = mesh.GetTriangle(i);

         for(size_t j = 0; j < qp.order; j++)
         {
            int idx = (i * qp.order + j) * 3;

            coords[idx + 0] = tr.XFromST(qp.x[j], qp.y[j]);
            coords[idx + 1] = tr.YFromST(qp.x[j], qp.y[j]);
            coords[idx + 2] = tr.ZFromST(qp.x[j], qp.y[j]);
         }

         areas[i] = tr.Area();
      }

      vector<double> weights(qp.order);

      for(size_t i = 0; i < qp.order; i++)
      {
         weights[i] = qp.w[i];
      }

      result.resize(points.size());

      calcIntegralOverArray(f,
                            points_double.data(),
                            coords.data(),
                            weights.data(),
                            areas.data(),
                            points.size(),
                            mesh.TriangleCount(),
                            qp.order,
                            result.data());
   }

   void calcIntegralOverArray(double (*f)(double, double, double, double),
                              const double* points,
                              const double* coords,
                              const double* weights,
                              const double* areas,
                              const int pointsCount,
                              const int trianglesCount,
                              const int quadratureOrder,
                              double* result)
   {
      for(size_t i = 0; i < pointsCount; i++)
      {
         double total_sum = 0;

         for(size_t j = 0; j < trianglesCount; j++)
         {
            double triangle_sum = 0;

            for(size_t k = 0; k < quadratureOrder; k++)
            {
               int idx = (j * quadratureOrder + k) * 3;
               triangle_sum += weights[k] * f(coords[idx + 0],
                                           coords[idx + 1],
                                           coords[idx + 2],
                                           points[i]);
            }

            total_sum += triangle_sum * areas[j];
         }

         result[i] = total_sum;
      }
   }
}