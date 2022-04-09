#include "quadrature.hpp"

Quadrature::Quadrature()
{
}

Quadrature::Quadrature(const Vector3& coordinates, real volume, real weight):
   coordinates(coordinates), weight(volume * weight)
{

}
