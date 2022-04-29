#include "testing_helpers.hpp"

namespace test
{

   Torus createTorus()
   {
      const double torusRadius = 2;
      const double torusSectionWidth = 0.2;
      
      //return Torus(torusRadius, torusSectionWidth, 80, 16, 16);
      return Torus(torusRadius, torusSectionWidth, 20, 4, 4);
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
      Vector3 step = (end - begin) / (n - 1);
      std::vector<Vector3> res(n);

      for(size_t i = 0; i < n; i++)
      {
         res[i] = begin + step * i;
      }

      return res;
   }
}
