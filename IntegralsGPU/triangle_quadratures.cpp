#include "triangle_quadratures.h"
#include <vector>

namespace triangle_quadratures
{
   float calcIntegralOverTriangle(float (*f)(Vector3),
                                   const Triangle& tr,
                                   const BasisQuadratures& qp)
   {
      float results = 0;

      for (size_t i = 0; i < qp.order; i++)
      {
         Vector3 v = tr.PointFromST(qp.x[i], qp.y[i]);
         results += qp.w[i] * f(v);
      }

      return results * tr.Area();
   }

   float calcIntegralOverMesh(float(*f)(Vector3),
                               const Mesh& mesh,
                               const BasisQuadratures& qp)
   {
      float integral_sum = 0;

      for (size_t t = 0; t < mesh.TrianglesCount(); t++)
      {
         integral_sum += calcIntegralOverTriangle(f, mesh.GetTriangle(t), qp);
      }

     return integral_sum;
   }

   float one(Vector3 v)
   {
      return 1;
   }

   float calcSurfaceArea(Mesh& mesh,
                          const BasisQuadratures& qp)
   {
      return calcIntegralOverMesh(one, mesh, qp);
   }

   /*void calcIntegralOverMesh(float(*f)(Vector3),
                             const Mesh& mesh, 
                             const BasisQuadratures& qp,
                             const vector<Vector3>& points,
                             std::vector<float>& results)
   {
      vector<float> points_float(points.size() * 3);

      for(size_t i = 0; i < points.size(); i++)
      {
         points_float[i * 3 + 0] = points[i].x;
         points_float[i * 3 + 1] = points[i].y;
         points_float[i * 3 + 2] = points[i].z;
      }

      vector<float> coords(mesh.TrianglesCount() * qp.order * 3);
      vector<float> areas(mesh.TrianglesCount());

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

      vector<float> weights(qp.order);

      for(size_t i = 0; i < qp.order; i++)
      {
         weights[i] = qp.w[i];
      }

      vector<float> normals(mesh.TrianglesCount() * 3);

      results.resize(points.size());

      calcIntegralOverArray(f,
                            points_float.data(),
                            coords.data(),
                            weights.data(),
                            areas.data(),
                            normals.data(),
                            points.size(),
                            mesh.TrianglesCount(),
                            qp.order,
                            results.data());
   }*/

//   void calcIntegralOverArray(float (*f)(float*),
//                              const float* points,
//                              const float* coords,
//                              const float* weights,
//                              const float* normals,
//                              const float* areas,
//                              const int pointsCount,
//                              const int trianglesCount,
//                              const int quadratureOrder,
//                              float* results)
//   {
//      float* params = new float[5];
//
//      for(size_t i = 0; i < pointsCount; i++)
//      {
//         float total_sum = 0;
//
//         for(size_t j = 0; j < trianglesCount; j++)
//         {
//            float triangle_sum = 0;
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
//         results[i] = total_sum;
//      }
//
//      delete[] params;
//   }
}