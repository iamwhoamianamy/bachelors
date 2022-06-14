#pragma once
#include "quadrature.hpp"

struct BEMQuadrature : public Quadrature
{
   Vector3 B;
   Vector3 normal;

   BEMQuadrature() = default;
   BEMQuadrature(
      const Vector3& coordinates,
      real area,
      real weight,
      const Vector3& B,
      const Vector3& normal);
};