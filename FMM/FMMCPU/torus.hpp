#pragma once
#include <vector>
#include "vector3.cuh"
#include "tetrahedron.hpp"
#include "hexahedron.hpp"

class Torus
{
private:
   const double radius;
   const double sectionWidth;
   const int onLengthStepCount;
   const int onWidthStepCount;
   const int onHeightStepCount;
public:
   Torus(const double radius,
         const double sectionWidth,
         const int onLengthStepCount,
         const int onWidthStepCount,
         const int onHeightStepCount);
   double innerRadius() const;
   double outerRadius() const;
   double bottom() const;
   double top() const;
   double stepAngle() const;
   std::vector<Tetrahedron> tetrahedra;

private:
   void buildHexahedra(std::vector<Hexahedron>& hexahedra);
   void buildTetrahedra(std::vector<Hexahedron>& hexahedra);
};

