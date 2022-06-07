#include "bem_quadrature.hpp"

BEMQuadrature::BEMQuadrature(
   const Vector3& coordinates,
   real volume,
   real weight,
   const Vector3& B):
   Quadrature(coordinates, volume, weight), B(B)
{

}
