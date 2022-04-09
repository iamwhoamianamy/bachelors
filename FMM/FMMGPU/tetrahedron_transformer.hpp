#pragma once
#include "vector3.cuh"
#include "tetrahedron.hpp"

class TetrahedronTransformer
{
public:
   static Vector3 pointFromBasisQuadrature(const Tetrahedron& tetr, const Vector3 basisQuadrature);
};

