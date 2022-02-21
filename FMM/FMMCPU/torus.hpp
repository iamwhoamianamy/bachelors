#pragma once
#include <vector>
#include "vector3.hpp"
#include "hexahedron.hpp"

class Torus
{
public:
   Torus(const double radius,
         const double sectionWidth,
         const int onLengthStepCount,
         const int onWidthStepCount,
         const int onHeightStepCount);

   std::vector<Hexahedron> hexahedrons;
};

