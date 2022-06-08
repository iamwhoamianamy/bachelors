#include "bem_quadrature.hpp"

BEMQuadrature::BEMQuadrature(
   const Vector3& coordinates,
   real square,
   real weight,
   const Vector3& B):
   Quadrature(coordinates, square, weight), B(B)
{

}
