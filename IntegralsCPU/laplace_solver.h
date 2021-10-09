#pragma once
#include <vector>
#include "triangle_quadratures.h"

using namespace triangle_quadratures;

class LaplaceSolver
{
public:
   double U(Vector3 point);
   Vector3 GradU(Vector3 point);

   double LaplaceIntegral1(Vector3 v,
                           Vector3 point,
                           Vector3 normal);

   double LaplaceIntegral2(Vector3 v,
                           Vector3 point,
                           Vector3 normal);

   void calcIntegralOverMesh(const Mesh& mesh,
                             const QuadPoints& qp,
                             const vector<Vector3>& points,
                             std::vector<double>& result)
};