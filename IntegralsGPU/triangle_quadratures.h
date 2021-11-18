#pragma once
#include <vector>
#include "triangle.h"
#include "exeptions.h"
#include "basis_quadratures.h"
#include "vector3.cuh"
#include "mesh.h"

namespace triangle_quadratures
{
   float calcIntegralOverTriangle(float (*f)(Vector3),
                                   const Triangle& tr,
                                   const BasisQuadratures& qp);

   float calcIntegralOverMesh(float (*f)(Vector3),
                               const Mesh& mesh,
                               const BasisQuadratures& qp);

   float calcSurfaceArea(Mesh& mesh,
                          const BasisQuadratures& qp);

   //void calcIntegralOverArray(float (*f)(float*),
   //                           const float* points,
   //                           const float* coords,
   //                           const float* weights,
   //                           const float* areas,
   //                           const float* normals,
   //                           const int pointsCount,
   //                           const int trianglesCount,
   //                           const int quadratureOrder,
   //                           float* results);
}