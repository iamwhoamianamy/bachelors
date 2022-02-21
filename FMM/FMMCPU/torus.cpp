#include <cmath>
#include "torus.hpp"
#include "hexahedron.hpp"

Torus::Torus(const double radius,
             const double sectionWidth,
             const int onLengthStepCount,
             const int onWidthStepCount,
             const int onHeightStepCount)
{
   const double innerRadius = radius - sectionWidth / 2.0;
   const double outerRadius = radius + sectionWidth / 2.0;
   const double bottom = -1.0 * sectionWidth / 2.0;
   const double top = sectionWidth / 2.0;
   const double stepAngle = 360.0 / onLengthStepCount * 3.14159265359 / 180.0;
   const double hexahedronHeight = sectionWidth / onHeightStepCount;

   auto getInnerRadiusPoint = [&](int n)
   {
      return Vector3(innerRadius * cos(stepAngle * n), innerRadius * sin(stepAngle * n), 0);
   };

   auto getPointOnWidthSection = [&](Vector3 innerRadiusPoint, int n)
   {
      double onWidthStepLength = sectionWidth / onWidthStepCount;
      return innerRadiusPoint + innerRadiusPoint.Normalized() * onWidthStepLength * n;
   };

   hexahedrons = std::vector<Hexahedron>(onLengthStepCount * onWidthStepCount * onHeightStepCount);
   size_t hexahedronsCount = 0;

   for(size_t lengthCounter = 0; lengthCounter < onLengthStepCount; lengthCounter++)
   {
      Vector3 innerRadiusPoint = getInnerRadiusPoint(lengthCounter);
      Vector3 nextInnerRadiusPoint = getInnerRadiusPoint(lengthCounter + 1);

      for(size_t widthCounter = 0; widthCounter < onWidthStepCount; widthCounter++)
      {
         Vector3 rightOnWidthSectionPoint = getPointOnWidthSection(innerRadiusPoint, widthCounter);
         Vector3 nextRightOnWidthSectionPoint = getPointOnWidthSection(innerRadiusPoint, widthCounter + 1);

         Vector3 leftOnWidthSectionPoint = getPointOnWidthSection(nextInnerRadiusPoint, widthCounter);
         Vector3 nextLeftOnWidthSectionPoint = getPointOnWidthSection(nextInnerRadiusPoint, widthCounter + 1);

         for(size_t heightCounter = 0; heightCounter < onHeightStepCount; heightCounter++)
         {
            double hexahedronBottom = bottom + hexahedronHeight * heightCounter;
            double hexahedronTop = bottom + hexahedronHeight * (heightCounter + 1);

            std::vector<Vector3> points(8);

            points[0] = Vector3(rightOnWidthSectionPoint.x, rightOnWidthSectionPoint.y, hexahedronBottom);
            points[1] = Vector3(nextRightOnWidthSectionPoint.x, nextRightOnWidthSectionPoint.y, hexahedronBottom);
            points[2] = Vector3(leftOnWidthSectionPoint.x, leftOnWidthSectionPoint.y, hexahedronBottom);
            points[3] = Vector3(nextLeftOnWidthSectionPoint.x, nextLeftOnWidthSectionPoint.y, hexahedronBottom);

            points[4] = Vector3(rightOnWidthSectionPoint.x, rightOnWidthSectionPoint.y, hexahedronTop);
            points[5] = Vector3(nextRightOnWidthSectionPoint.x, nextRightOnWidthSectionPoint.y, hexahedronTop);
            points[6] = Vector3(leftOnWidthSectionPoint.x, leftOnWidthSectionPoint.y, hexahedronTop);
            points[7] = Vector3(nextLeftOnWidthSectionPoint.x, nextLeftOnWidthSectionPoint.y, hexahedronTop);

            hexahedrons[hexahedronsCount++] = Hexahedron(points);
         }
      }
   }
   
}



