#include "cylinder.hpp"

#include "math.hpp"

Cylinder::Cylinder(
   real height,
   real radius,
   size_t widthSegmentCount,
   size_t heightSegmentCount) :
   _height(height),
   _radius(radius),
   _widthSegmentCount(widthSegmentCount),
   _heightSegmentCount(heightSegmentCount)
{
   _sideTriangles.reserve(widthSegmentCount * heightSegmentCount * 2);

   real width = 2 * math::PI;
   real segmentWidth = width / widthSegmentCount;
   real segmentHeight = height / heightSegmentCount;

   Vector3 widthSegmentStep(segmentWidth, 0, 0);
   Vector3 heightSegmentStep(0, 0, segmentHeight);

   for(size_t h = 0; h < heightSegmentCount; h++)
   {
      Vector3 currentHeight(0, 0, segmentHeight * h);

      for(size_t w = 0; w < widthSegmentCount; w++)
      {
         Vector3 currentWidth(segmentWidth * w, 0, 0);
         Vector3 triangleBase = currentWidth + currentHeight;

         _sideTriangles.emplace_back(
            triangleBase,
            triangleBase + widthSegmentStep,
            triangleBase + heightSegmentStep);

         _sideTriangles.emplace_back(
            triangleBase + widthSegmentStep + heightSegmentStep,
            triangleBase + widthSegmentStep,
            triangleBase + heightSegmentStep);
      }
   }
}

std::vector<Triangle>& Cylinder::sideTriangles()
{
   return _sideTriangles;
}

const std::vector<Triangle>& Cylinder::sideTriangles() const
{
   return _sideTriangles;
}
