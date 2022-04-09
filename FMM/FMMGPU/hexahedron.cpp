#include "hexahedron.hpp"

Hexahedron::Hexahedron()
{
   points.resize(8);
}

Hexahedron::Hexahedron(std::vector<Vector3> points)
   :points(points)
{

}