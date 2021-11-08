#pragma once
#include <vector>
#include "vector3.cuh"
#include "mesh.h"
#include "quad_points.h"

using namespace triangle_quadratures;
using namespace std;

namespace laplace_solver
{
   double u(Vector3 point);
   Vector3 gradU(Vector3 point);

   double laplaceIntegral1(Vector3 v,
                           Vector3 point,
                           Vector3 normal);

   double laplaceIntegral2(Vector3 v,
                           Vector3 point,
                           Vector3 normal);

   void calcIntegralOverMesh(const Mesh& mesh,
                             const QuadPoints& qp,
                             const vector<Vector3>& points,
                             std::vector<double>& result);
};