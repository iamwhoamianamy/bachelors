#pragma once
#include <vector>
#include "triangle.h"
#include "exeptions.h"
#include "quad_points.h"
#include "point.h"
#include "mesh.h"

namespace triangle_quadratures
{
   double calcIntegralOverTriangle(double (*f)(Point),
                                   const Triangle& tr,
                                   const QuadPoints& qp);

   void calcIntegralOverMesh(double (*f)(Point),
                             const Mesh& mesh,
                             const QuadPoints& qp,
                             const vector<Point>& points,
                             std::vector<double>& result);
}