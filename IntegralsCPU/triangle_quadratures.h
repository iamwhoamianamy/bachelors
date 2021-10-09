#pragma once
#include <vector>
#include "triangle.h"
#include "exeptions.h"
#include "quad_points.h"
#include "vector3.h"
#include "mesh.h"

namespace triangle_quadratures
{
   double calcIntegralOverTriangle(double (*f)(Vector3 v),
                                   const Triangle& tr,
                                   const QuadPoints& qp);

   void calcIntegralOverMesh(double (*f)(Vector3 v),
                             const Mesh& mesh,
                             const QuadPoints& qp,
                             const vector<Vector3>& points,
                             std::vector<double>& result);

   //void calcIntegralOverArray(double (*f)(double*),
   //                           const double* points,
   //                           const double* coords,
   //                           const double* weights,
   //                           const double* areas,
   //                           const double* normals,
   //                           const int pointsCount,
   //                           const int trianglesCount,
   //                           const int quadratureOrder,
   //                           double* result);
}