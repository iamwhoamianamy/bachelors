#include "cylinder.hpp"

#include "math.hpp"

Cylinder::Cylinder(
   real bottom,
   real top,
   real radius,
   size_t widthSegmentCount,
   size_t heightSegmentCount,
   size_t depthSegmentCount)
{
   _sideTriangles.reserve(widthSegmentCount * heightSegmentCount * 2);

   real width = 2 * math::PI;
   real segmentWidth = width / widthSegmentCount;
   real segmentHeight = (top - bottom) / heightSegmentCount;

   Vector3 widthSegmentStep(segmentWidth, 0, 0);
   Vector3 heightSegmentStep(0, 0, segmentHeight);

   for(size_t h = 0; h < heightSegmentCount; h++)
   {
      Vector3 currentHeight(0, 0, segmentHeight * h + bottom);

      for(size_t w = 0; w < widthSegmentCount; w++)
      {
         Vector3 currentWidth(segmentWidth * w, 0, 0);
         Vector3 triangleBase = currentWidth + currentHeight + Vector3(0, radius, 0);

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

   real segmentDepth = radius / depthSegmentCount;
   Vector3 depthSegmentStep(0, segmentDepth, 0);

   _topTriangles.reserve(depthSegmentCount * widthSegmentCount * 2);

   for(size_t d = 0; d < depthSegmentCount; d++)
   {
      Vector3 currentDepth(0, segmentDepth * d, 0);

      for(size_t w = 0; w < widthSegmentCount; w++)
      {
         Vector3 currentWidth(segmentWidth * w, 0, 0);
         Vector3 triangleBase = currentWidth + currentDepth + Vector3(0, 0, top);

         if(d == 0)
         {
            _topTriangles.emplace_back(
               triangleBase,
               triangleBase + depthSegmentStep,
               triangleBase + widthSegmentStep + depthSegmentStep);
         }
         else
         {
            _topTriangles.emplace_back(
               triangleBase,
               triangleBase + widthSegmentStep,
               triangleBase + depthSegmentStep);

            _topTriangles.emplace_back(
               triangleBase + widthSegmentStep + depthSegmentStep,
               triangleBase + widthSegmentStep,
               triangleBase + depthSegmentStep);
         }
      }
   }

   _bottomTriangles = _topTriangles;

   for (auto &bottomTriangle : _bottomTriangles)
   {
      bottomTriangle.a().z = bottom;
      bottomTriangle.b().z = bottom;
      bottomTriangle.c().z = bottom;
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

std::vector<Triangle>& Cylinder::topTriangles()
{
   return _topTriangles;
}

const std::vector<Triangle>& Cylinder::topTriangles() const
{
   return _topTriangles;
}

std::vector<Triangle>& Cylinder::bottomTriangles()
{
   return _bottomTriangles;
}

const std::vector<Triangle>& Cylinder::bottomTriangles() const
{
   return _bottomTriangles;
}
