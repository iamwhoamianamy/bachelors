#pragma once
#include <vector>
#include "vector3.cuh"
#include "mesh.h"
#include "basis_quadratures.h"

using namespace triangle_quadratures;
using namespace std;

namespace laplace_solver
{
   float u(Vector3 point);
   Vector3 gradU(Vector3 point);

   float laplaceIntegral1(Vector3 v,
                           Vector3 point,
                           Vector3 normal);

   float laplaceIntegral2(Vector3 v,
                           Vector3 point,
                           Vector3 normal);

   void calcIntegralOverMesh(const Mesh& mesh,
                             const BasisQuadratures& qp,
                             const vector<Vector3>& points,
                             std::vector<float>& results);
};