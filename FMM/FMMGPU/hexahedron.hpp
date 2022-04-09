#pragma once
#include "vector3.cuh"
#include <vector>

class Hexahedron
{
public:
   std::vector<Vector3> points;
   Hexahedron();
   Hexahedron(std::vector<Vector3> points);
};