#include "triangle_quadratures.h"
#include <vector>

namespace triangle_quadratures
{
   double calcIntegralOverTriangle(double (*f)(Vector3),
                                   const Triangle& tr,
                                   const QuadPoints& qp)
   {
      double result = 0;

      for (size_t i = 0; i < qp.order; i++)
      {
         Vector3 v = tr.PointFromST(qp.x[i], qp.y[i]);
         result += qp.w[i] * f(v);
      }

      return result * tr.Area();
   }

   double calcIntegralOverMesh(double(*f)(Vector3),
                               const Mesh& mesh,
                               const QuadPoints& qp)
   {
      double integral_sum = 0;

      for (size_t t = 0; t < mesh.TrianglesCount(); t++)
      {
         integral_sum += calcIntegralOverTriangle(f, mesh.GetTriangle(t), qp);
      }

     return integral_sum;
   }

   double one(Vector3 v)
   {
      return 1;
   }

   double calcSurfaceArea(Mesh& mesh,
                          const QuadPoints& qp)
   {
      return calcIntegralOverMesh(one, mesh, qp);
   }

   /*void calcIntegralOverMesh(double(*f)(Vector3),
                             const Mesh& mesh, 
                             const QuadPoints& qp,
                             const vector<Vector3>& points,
                             std::vector<double>& result)
   {
      vector<double> points_double(points.size() * 3);

      for(size_t i = 0; i < points.size(); i++)
      {
         points_double[i * 3 + 0] = points[i].x;
         points_double[i * 3 + 1] = points[i].y;
         points_double[i * 3 + 2] = points[i].z;
      }

      vector<double> coords(mesh.TrianglesCount() * qp.order * 3);
      vector<double> areas(mesh.TrianglesCount());

      for(size_t i = 0; i < mesh.TrianglesCount(); i++)
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

      vector<double> normals(mesh.TrianglesCount() * 3);

      result.resize(points.size());

      calcIntegralOverArray(f,
                            points_double.data(),
                            coords.data(),
                            weights.data(),
                            areas.data(),
                            normals.data(),
                            points.size(),
                            mesh.TrianglesCount(),
                            qp.order,
                            result.data());
   }*/

//   void calcIntegralOverArray(double (*f)(double*),
//                              const double* points,
//                              const double* coords,
//                              const double* weights,
//                              const double* normals,
//                              const double* areas,
//                              const int pointsCount,
//                              const int trianglesCount,
//                              const int quadratureOrder,
//                              double* result)
//   {
//      double* params = new double[5];
//
//      for(size_t i = 0; i < pointsCount; i++)
//      {
//         double total_sum = 0;
//
//         for(size_t j = 0; j < trianglesCount; j++)
//         {
//            double triangle_sum = 0;
//
//            for(size_t k = 0; k < quadratureOrder; k++)
//            {
//               int idx = (j * quadratureOrder + k) * 3;
//
//               params[0] = coords[idx + 0];
//               params[0] =    coords[idx + 1],
//               params[0] =    coords[idx + 2],
//               params[0] =    points[i];
//               params[0] =    points[i];
//
//               triangle_sum += weights[k] * f(params.data());
//            }
//
//            total_sum += triangle_sum * areas[j];
//         }
//
//         result[i] = total_sum;
//      }
//
//      delete[] params;
//   }
}