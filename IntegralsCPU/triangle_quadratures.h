#pragma once
#include <vector>
#include "triangle.h"
#include "exeptions.h"
#include "quad_points.h"
#include "point.h"
#include "mesh.h"

namespace triangle_quadratures
{
   double calcIntegralOverTriangle(double (*f)(double, double, double, double),
                                   const Triangle& tr,
                                   const QuadPoints& qp);

   void calcIntegralOverMesh(double (*f)(double, double, double, double),
                             const Mesh& mesh,
                             const QuadPoints& qp,
                             const vector<Point>& points,
                             std::vector<double>& result);

   void calcIntegralOverArray(double (*f)(double, double, double, double),
                              const double* points,
                              const double* coords,
                              const double* weights,
                              const double* areas,
                              const int pointsCount,
                              const int trianglesCount,
                              const int quadratureOrder,
                              double* result);
}