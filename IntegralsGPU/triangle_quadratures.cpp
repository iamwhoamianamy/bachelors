#include "triangle_quadratures.h"
#include <vector>

namespace triangle_quadratures
{
   real calcIntegralOverTriangle(real (*f)(Vector3),
                                   const Triangle& tr,
                                   const BasisQuadratures& qp)
   {
      real results = 0;

      for (size_t i = 0; i < qp.order; i++)
      {
         Vector3 v = tr.PointFromST(qp.x[i], qp.y[i]);
         results += qp.w[i] * f(v);
      }

      return results * tr.Area();
   }

   real calcIntegralOverMesh(real(*f)(Vector3),
                               const Mesh& mesh,
                               const BasisQuadratures& qp)
   {
      real integral_sum = 0;

      for (size_t t = 0; t < mesh.TrianglesCount(); t++)
      {
         integral_sum += calcIntegralOverTriangle(f, mesh.GetTriangle(t), qp);
      }

     return integral_sum;
   }

   real one(Vector3 v)
   {
      return 1;
   }

   real calcSurfaceArea(Mesh& mesh,
                          const BasisQuadratures& qp)
   {
      return calcIntegralOverMesh(one, mesh, qp);
   }

   /*void calcIntegralOverMesh(real(*f)(Vector3),
                             const Mesh& mesh, 
                             const BasisQuadratures& qp,
                             const vector<Vector3>& points,
                             std::vector<real>& results)
   {
      vector<real> points_real(points.size() * 3);

      for(size_t i = 0; i < points.size(); i++)
      {
         points_real[i * 3 + 0] = points[i].x;
         points_real[i * 3 + 1] = points[i].y;
         points_real[i * 3 + 2] = points[i].z;
      }

      vector<real> coords(mesh.TrianglesCount() * qp.order * 3);
      vector<real> areas(mesh.TrianglesCount());

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

      vector<real> weights(qp.order);

      for(size_t i = 0; i < qp.order; i++)
      {
         weights[i] = qp.w[i];
      }

      vector<real> normals(mesh.TrianglesCount() * 3);

      results.resize(points.size());

      calcIntegralOverArray(f,
                            points_real.data(),
                            coords.data(),
                            weights.data(),
                            areas.data(),
                            normals.data(),
                            points.size(),
                            mesh.TrianglesCount(),
                            qp.order,
                            results.data());
   }*/

//   void calcIntegralOverArray(real (*f)(real*),
//                              const real* points,
//                              const real* coords,
//                              const real* weights,
//                              const real* normals,
//                              const real* areas,
//                              const int pointsCount,
//                              const int trianglesCount,
//                              const int quadratureOrder,
//                              real* results)
//   {
//      real* params = new real[5];
//
//      for(size_t i = 0; i < pointsCount; i++)
//      {
//         real total_sum = 0;
//
//         for(size_t j = 0; j < trianglesCount; j++)
//         {
//            real triangle_sum = 0;
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