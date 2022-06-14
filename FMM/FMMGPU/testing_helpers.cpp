#include <iostream>
#include <fstream>

#include "testing_helpers.hpp"
#include "math.hpp"
#include "exeptions.hpp"
#include "integration.hpp"

namespace test
{
   Torus createTorus()
   {
      constexpr real torusRadius = 2;
      constexpr real torusSectionWidth = 0.2;
      
      //return { torusRadius, torusSectionWidth, 80, 8, 8 };
      //return { torusRadius, torusSectionWidth, 40, 8, 8 };
      return { torusRadius, torusSectionWidth, 20, 4, 4 };
   }

   Cylinder createCylinder()
   {
      real cylinderRadius = 1.0;
      real cylinderBottom = -1.5;
      real cylinderTop = 1.5;

      real subdivisionLevel = 1;

      size_t widthSegmentCount = 128 * subdivisionLevel;
      size_t heightSegmentCount = 99 * subdivisionLevel;
      size_t depthSegmentCount = 50 * subdivisionLevel;

      return {
         cylinderBottom,
         cylinderTop,
         cylinderRadius,
         widthSegmentCount,
         heightSegmentCount,
         depthSegmentCount };
   }
   
   std::vector<ReferenceCylinderData> readCylinderData(const std::string& filename)
   {
      std::vector<ReferenceCylinderData> result;

      std::ifstream fin(filename);
      std::string _;
      std::getline(fin, _);
      std::getline(fin, _);
      size_t pointId;

      while(fin >> pointId)
      {
         real px, py, pz;
         real bx, by, bz, bl;

         fin >> px >> py >> pz >> bx >> by >> bz >> bl;

         result.emplace_back(pointId, Vector3(px, py, pz), Vector3(bx, by, bz), bl);
      }

      return result;
   }

   std::vector<BEMQuadrature> quadraturesFromCylinder()
   {
      Cylinder cylinder = createCylinder();
      BasisQuadratures bq = readTriangleBasisQuadratures();

      auto externalCylinderSides = readCylinderData("cylinder/¬нешний÷илиндр.0");
      auto externalCylinderTop = readCylinderData("cylinder/¬нешний÷илиндр¬ерх.0");
      auto externalCylinderBottom = readCylinderData("cylinder/¬нешний÷илиндрЌиз.0");

      auto BEMQuadraturesSide = math::calcBEMquadraturesFromTriangles(
         cylinder.sideTriangles(), bq, externalCylinderSides, 0);

      auto BEMQuadraturesTop = math::calcBEMquadraturesFromTriangles(
         cylinder.topTriangles(), bq, externalCylinderTop, 1);

      auto BEMQuadraturesBottom = math::calcBEMquadraturesFromTriangles(
         cylinder.bottomTriangles(), bq, externalCylinderBottom, -1);

      std::vector<BEMQuadrature> quadratures;
      quadratures.reserve(
         BEMQuadraturesSide.size() +
         BEMQuadraturesTop.size() +
         BEMQuadraturesTop.size());

      quadratures.insert(quadratures.end(), BEMQuadraturesSide.begin(), BEMQuadraturesSide.end());
      quadratures.insert(quadratures.end(), BEMQuadraturesTop.begin(), BEMQuadraturesTop.end());
      quadratures.insert(quadratures.end(), BEMQuadraturesBottom.begin(), BEMQuadraturesBottom.end());

      return quadratures;
   }

   double getTime(void (*f)())
   {
      auto start = std::chrono::steady_clock::now();
      f();
      auto stop = std::chrono::steady_clock::now();
      return getTime(start, stop);
   }

   double getTime(const std::chrono::steady_clock::time_point& start,
                  const std::chrono::steady_clock::time_point& stop)
   {
      return std::chrono::duration_cast<std::chrono::microseconds>
         (stop - start).count() * 1e-6;
   }

   BasisQuadratures readTetrahedronBasisQuadratures()
   {
      BasisQuadratures bq;
      std::string bqDir = "E:/ћнЄ/Ѕиба/bachelors/FMM/FMMGPU/basis_quadratures/";

      try
      {
         bq.initFromTXT(bqDir + "keast 7 x.txt", bqDir + "keast 7 w.txt");
         //std::cout << "Quadratures were successfully read. Order = " << bq._order() << std::endl;
      }
      catch(Exeption ex)
      {
         std::cout << ex;
      }

      return bq;
   }

   BasisQuadratures readTriangleBasisQuadratures()
   {
      BasisQuadratures bq;
      std::string path = "basis_quadratures/";

      try
      {
         bq.initFromTXT(path + "gauss7_xy.txt", path + "gauss7_w.txt");
      }
      catch(Exeption ex)
      {
         std::cout << ex.message;
      }

      return bq;
   }

   std::vector<Vector3> createPoints(const Vector3& begin, const Vector3& end, int n)
   {
      if(n == 1)
      {
         return { (end + begin) / 2 };
      }
      else
      {
         Vector3 step = (end - begin) / (n - 1);
         std::vector<Vector3> res(n);

         for(size_t i = 0; i < n; i++)
         {
            res[i] = begin + step * i;
         }

         return res;
      }
   }



   std::vector<Vector3> createRandomPoints(const Box& box, int n)
   {
      std::vector<Vector3> res(n);
      std::srand(std::time(0));

      for (int i = 0; i < n; ++i)
      {
         real x = math::randBetween(
            box.center().x - box.halfDimensions().x,
            box.center().x + box.halfDimensions().x);

         real y = math::randBetween(
            box.center().y - box.halfDimensions().y,
            box.center().y + box.halfDimensions().y);

         real z = math::randBetween(
            box.center().z - box.halfDimensions().z,
            box.center().z + box.halfDimensions().z);

         res[i] = { x, y, z };
      }

      return res;
   }
   

   void printSeparateLine(std::ostream& os, size_t count)
   {
      for (int i = 0; i < count; ++i)
      {
         os << "-";
      }

      os << std::endl;
   }
}
