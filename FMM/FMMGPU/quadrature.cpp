#include "quadrature.hpp"

Quadrature::Quadrature()
{
}

Quadrature::Quadrature(const Vector3& coordinates, real volume, real weight):
   Vector3(coordinates), weight(volume * weight)
{

}
