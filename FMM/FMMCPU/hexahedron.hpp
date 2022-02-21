#pragma once
#include "vector3.hpp"
#include <vector>

class Hexahedron
{
private:
   std::vector<Vector3> points;
public:
   Hexahedron();
   Hexahedron(std::vector<Vector3> points);
};