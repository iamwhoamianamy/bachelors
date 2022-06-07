#pragma once
#include "real.hpp"
#include "triangle.hpp"

class Cylinder
{
private:
   real _height;
   real _radius;
   size_t _widthSegmentCount;
   size_t _heightSegmentCount;
   std::vector<Triangle> _sideTriangles;

public:
   Cylinder(
      real height,
      real radius,
      size_t widthSegmentCount,
      size_t heightSegmentCount);

   std::vector<Triangle>& sideTriangles();
   const std::vector<Triangle>& sideTriangles() const;

private:
};
