#pragma once
#include "quadrature.hpp"

struct BEMQuadrature : public Quadrature
{
   Vector3 B;

   BEMQuadrature() = default;
   BEMQuadrature(
      const Vector3& coordinates,
      real square,
      real weight,
      const Vector3& B);
};