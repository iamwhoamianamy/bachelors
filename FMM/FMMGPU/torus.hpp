#pragma once
#include <vector>
#include "vector3.cuh"
#include "tetrahedron.hpp"
#include "hexahedron.hpp"
#include "real.hpp"

class Torus
{
private:
   const real _radius;
   const real _sectionWidth;
   const int _onLengthStepCount;
   const int _onWidthStepCount;
   const int _onHeightStepCount;
public:
   Torus(const real radius,
         const real sectionWidth,
         const int onLengthStepCount,
         const int onWidthStepCount,
         const int onHeightStepCount);
   real innerRadius() const;
   real outerRadius() const;
   real bottom() const;
   real top() const;
   real stepAngle() const;
   std::vector<Tetrahedron> tetrahedra;

private:
   void buildHexahedra(std::vector<Hexahedron>& hexahedra);
   void buildTetrahedra(std::vector<Hexahedron>& hexahedra);
};

