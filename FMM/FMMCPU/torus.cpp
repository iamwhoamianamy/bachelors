#include <cmath>
#include "torus.hpp"
#include "hexahedron.hpp"
#include <fstream>

void saveHexahedraAsOBJ(std::vector<Hexahedron>& hexahedra, std::string filename);
void saveTetrahedraAsOBJ(std::vector<Tetrahedron>& tetrahedra, std::string filename);

Torus::Torus(const double radius,
             const double sectionWidth,
             const int onLengthStepCount,
             const int onWidthStepCount,
             const int onHeightStepCount)
   : radius(radius), sectionWidth(sectionWidth),
   onLengthStepCount(onLengthStepCount),
   onWidthStepCount(onWidthStepCount),
   onHeightStepCount(onHeightStepCount)
{
   std::vector<Hexahedron> hexahedra(onLengthStepCount * onWidthStepCount * onHeightStepCount);
   buildHexahedra(hexahedra);
   buildTetrahedra(hexahedra);

   saveHexahedraAsOBJ(hexahedra, "testmesh/debsHex.obj");
   saveTetrahedraAsOBJ(tetrahedra, "testmesh/debsTet.obj");
}

double Torus::innerRadius() const
{
   return radius - sectionWidth / 2.0;
}

double Torus::outerRadius() const
{
   return radius + sectionWidth / 2.0;
}

double Torus::bottom() const
{
   return -1.0 * sectionWidth / 2.0;
}

double Torus::top() const
{
   return sectionWidth / 2.0;
}

double Torus::stepAngle() const
{
   return 360.0 / onLengthStepCount * 3.14159265359 / 180.0;
}

void Torus::buildHexahedra(std::vector<Hexahedron>& hexahedra)
{
   const double hexahedronHeight = sectionWidth / onHeightStepCount;

   auto getInnerRadiusPoint = [&](int n)
   {
      return Vector3(innerRadius() * cos(stepAngle() * n), innerRadius() * sin(stepAngle() * n), 0);
   };

   auto getPointOnWidthSection = [&](Vector3 innerRadiusPoint, int n)
   {
      double onWidthStepLength = sectionWidth / onWidthStepCount;
      return innerRadiusPoint + innerRadiusPoint.Normalized() * onWidthStepLength * n;
   };

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
            double hexahedronBottom = bottom() + hexahedronHeight * heightCounter;
            double hexahedronTop = bottom() + hexahedronHeight * (heightCounter + 1);

            std::vector<Vector3> points(8);

            points[0] = Vector3(rightOnWidthSectionPoint.x, rightOnWidthSectionPoint.y, hexahedronBottom);
            points[1] = Vector3(nextRightOnWidthSectionPoint.x, nextRightOnWidthSectionPoint.y, hexahedronBottom);
            points[2] = Vector3(leftOnWidthSectionPoint.x, leftOnWidthSectionPoint.y, hexahedronBottom);
            points[3] = Vector3(nextLeftOnWidthSectionPoint.x, nextLeftOnWidthSectionPoint.y, hexahedronBottom);

            points[4] = Vector3(rightOnWidthSectionPoint.x, rightOnWidthSectionPoint.y, hexahedronTop);
            points[5] = Vector3(nextRightOnWidthSectionPoint.x, nextRightOnWidthSectionPoint.y, hexahedronTop);
            points[6] = Vector3(leftOnWidthSectionPoint.x, leftOnWidthSectionPoint.y, hexahedronTop);
            points[7] = Vector3(nextLeftOnWidthSectionPoint.x, nextLeftOnWidthSectionPoint.y, hexahedronTop);

            hexahedra[hexahedronsCount++] = Hexahedron(points);
         }
      }
   }
}

void Torus::buildTetrahedra(std::vector<Hexahedron>& hexahedra)
{
   tetrahedra.reserve(hexahedra.size() * 5);

   for(auto& hex : hexahedra)
   {
      tetrahedra.emplace_back(Tetrahedron(hex.points[0], hex.points[1], hex.points[3], hex.points[5]));
      tetrahedra.emplace_back(Tetrahedron(hex.points[0], hex.points[3], hex.points[5], hex.points[6]));
      tetrahedra.emplace_back(Tetrahedron(hex.points[0], hex.points[2], hex.points[3], hex.points[6]));
      tetrahedra.emplace_back(Tetrahedron(hex.points[0], hex.points[4], hex.points[5], hex.points[6]));
      tetrahedra.emplace_back(Tetrahedron(hex.points[3], hex.points[5], hex.points[6], hex.points[7]));
   }
}

void saveHexahedraAsOBJ(std::vector<Hexahedron>& hexahedra, std::string filename)
{
   std::ofstream fout(filename);

   for(size_t i = 0; i < hexahedra.size(); i++)
   {
      for(auto& p : hexahedra[i].points)
      {
         fout << "v" << " " << p.x << " " << p.y << " " << p.z << std::endl;
      }
   }

   for(size_t i = 0; i < hexahedra.size(); i++)
   {
      int c = i * 8 + 1;

      fout << "f" << " " << c + 0 << " " << c + 1 << " " << c + 4 << std::endl; // right
      fout << "f" << " " << c + 1 << " " << c + 4 << " " << c + 5 << std::endl; // right

      fout << "f" << " " << c + 0 << " " << c + 2 << " " << c + 6 << std::endl; // inner
      fout << "f" << " " << c + 0 << " " << c + 4 << " " << c + 6 << std::endl; // inner

      fout << "f" << " " << c + 0 << " " << c + 1 << " " << c + 3 << std::endl; // bot
      fout << "f" << " " << c + 0 << " " << c + 2 << " " << c + 3 << std::endl; // bot

      fout << "f" << " " << c + 4 << " " << c + 5 << " " << c + 7 << std::endl; // top
      fout << "f" << " " << c + 4 << " " << c + 6 << " " << c + 7 << std::endl; // top

      fout << "f" << " " << c + 1 << " " << c + 3 << " " << c + 7 << std::endl; // outer
      fout << "f" << " " << c + 1 << " " << c + 5 << " " << c + 7 << std::endl; // outer

      fout << "f" << " " << c + 2 << " " << c + 6 << " " << c + 7 << std::endl; // left
      fout << "f" << " " << c + 2 << " " << c + 3 << " " << c + 7 << std::endl; // left
   }

   fout.close();
}

void saveTetrahedraAsOBJ(std::vector<Tetrahedron>& tetrahedra, std::string filename)
{
   std::ofstream fout(filename);

   for(size_t i = 0; i < tetrahedra.size(); i++)
   {
      for(auto& t : tetrahedra[i].points)
      {
         fout << "v" << " " << t.x << " " << t.y << " " << t.z << std::endl;
      }
   }

   for(size_t i = 0; i < tetrahedra.size(); i++)
   {
      int c = i * 4 + 1;

      fout << "f" << " " << c + 0 << " " << c + 1 << " " << c + 3 << std::endl;
      fout << "f" << " " << c + 0 << " " << c + 1 << " " << c + 2 << std::endl;
      fout << "f" << " " << c + 0 << " " << c + 2 << " " << c + 3 << std::endl;
      fout << "f" << " " << c + 1 << " " << c + 2 << " " << c + 3 << std::endl;
   }

   fout.close();
}


