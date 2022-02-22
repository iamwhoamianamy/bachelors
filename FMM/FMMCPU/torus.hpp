#pragma once
#include <vector>
#include "vector3.hpp"
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
   std::vector<Tetrahedron> tetrahedra;
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

private:
   void buildHexahedra(std::vector<Hexahedron>& hexahedra);
   void buildTetrahedra(std::vector<Hexahedron>& hexahedra);
};

