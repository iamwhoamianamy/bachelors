#pragma once
#include "real.h"
#include <vector>
#include "vector3.cuh"
#include "mesh.h"
#include "basis_quadratures.h"

using namespace triangle_quadratures;
using namespace std;

namespace laplace_solver
{
   real u(Vector3 point);
   Vector3 gradU(Vector3 point);

   real laplaceIntegral1(Vector3 v,
                         Vector3 point,
                         Vector3 normal);

   real laplaceIntegral2(Vector3 v,
                         Vector3 point,
                         Vector3 normal);

   void calcIntegralOverMesh(const Mesh& mesh,
                             const BasisQuadratures& qp,
                             const vector<Vector3>& points,
                             std::vector<real>& results);
};