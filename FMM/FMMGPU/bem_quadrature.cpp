#include "bem_quadrature.hpp"

BEMQuadrature::BEMQuadrature(
   const Vector3& coordinates,
   real area,
   real weight,
   const Vector3& B,
   const Vector3& normal):
   Quadrature(coordinates, area, weight), B(B), normal(normal)
{

}
