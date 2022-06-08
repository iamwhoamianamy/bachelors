#pragma once
#include "real.hpp"
#include "triangle.hpp"

class Cylinder
{
private:
   std::vector<Triangle> _sideTriangles;
   std::vector<Triangle> _topTriangles;
   std::vector<Triangle> _bottomTriangles;

public:
   Cylinder(
      real bottom,
      real top,
      real radius,
      size_t widthSegmentCount,
      size_t heightSegmentCount,
      size_t depthSegmentCount);

   std::vector<Triangle>& sideTriangles();
   const std::vector<Triangle>& sideTriangles() const;

   std::vector<Triangle>& topTriangles();
   const std::vector<Triangle>& topTriangles() const;

   std::vector<Triangle>& bottomTriangles();
   const std::vector<Triangle>& bottomTriangles() const;

private:
};
