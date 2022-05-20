#include <iostream>

#include "testing_helpers.hpp"
#include "math.hpp"
#include "exeptions.hpp"

namespace test
{
   Torus createTorus()
   {
      constexpr double torusRadius = 2;
      constexpr double torusSectionWidth = 0.2;
      
      //return { torusRadius, torusSectionWidth, 80, 8, 8 };
      //return { torusRadius, torusSectionWidth, 40, 8, 8 };
      return { torusRadius, torusSectionWidth, 20, 4, 4 };
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

   BasisQuadratures readBasisQuadratures()
   {
      BasisQuadratures bq;
      std::string bqDir = "E:/ÃÌ∏/¡Ë·‡/bachelors/FMM/FMMGPU/basis_quadratures/";
      //std::string bqDir = "basis_quadratures/";

      try
      {
         bq.InitFromTXT(bqDir + "keast 7 x.txt", bqDir + "keast 7 w.txt");
         //std::cout << "Quadratures were successfully read. Order = " << bq._order() << std::endl;
      }
      catch(Exeption ex)
      {
         std::cout << ex;
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
