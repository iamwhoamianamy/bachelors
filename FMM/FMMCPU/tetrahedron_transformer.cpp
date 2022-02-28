#include "tetrahedron_transformer.hpp"

Vector3 TetrahedronTransformer::pointFromBasisQuadrature(const Tetrahedron& tetr, const Vector3 basisQuadrature)
{
   return tetr.a() + 
      (tetr.b() - tetr.a()) * basisQuadrature.x +
      (tetr.c() - tetr.a()) * basisQuadrature.y + 
      (tetr.d() - tetr.a()) * basisQuadrature.z;
}
