#pragma once
#include "vector3.cuh"

struct Quadrature : public Vector3
{
   real weight;

   Quadrature();
   Quadrature(const Vector3& coordinates, real volume, real weight);
};