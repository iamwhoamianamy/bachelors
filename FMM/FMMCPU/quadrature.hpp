#pragma once
#include "vector3.hpp"

struct Quadrature
{
   Vector3 coordinates;
   real weight;

   Quadrature();
   Quadrature(const Vector3& coordinates, real volume, real weight);
};