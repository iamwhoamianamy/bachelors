#pragma once
#include "vector3.hpp"
#include "tetrahedron.hpp"

class TetrahedronTransformer
{
public:
   static Vector3 pointFromBasisQuadrature(const Tetrahedron& tetr, const Vector3 basisQuadrature);
};

