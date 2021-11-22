#pragma once
#include "real.h"
#include <vector>
#include "triangle.h"
#include "exeptions.h"
#include "basis_quadratures.h"
#include "vector3.cuh"
#include "mesh.h"

namespace triangle_quadratures
{
   real calcIntegralOverTriangle(real (*f)(Vector3),
                                 const Triangle& tr,
                                 const BasisQuadratures& qp);

   real calcIntegralOverMesh(real (*f)(Vector3),
                             const Mesh& mesh,
                             const BasisQuadratures& qp);

   real calcSurfaceArea(Mesh& mesh,
                        const BasisQuadratures& qp);

   //void calcIntegralOverArray(real (*f)(real*),
   //                           const real* points,
   //                           const real* coords,
   //                           const real* weights,
   //                           const real* areas,
   //                           const real* normals,
   //                           const int pointsCount,
   //                           const int trianglesCount,
   //                           const int quadratureOrder,
   //                           real* results);
}