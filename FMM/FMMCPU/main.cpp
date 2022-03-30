#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <chrono>

#include "vector3.hpp"
#include "tetrahedron.hpp"
#include "hexahedron.hpp"
#include "torus.hpp"
#include "basis_quadratures.hpp"
#include "exeptions.hpp"
#include "math.hpp"
#include "harmonics.hpp"
#include "multipole_solver.hpp"


using namespace math;
const real current = 5;

std::vector<std::pair<Vector3, Vector3>> readTelmaResults(const std::string& filename)
{
   std::vector<std::pair<Vector3, Vector3>> res;
   std::ifstream fin(filename);
   size_t pointsCount;
   std::string _;
   fin >> _ >> _ >> _;
   fin >> pointsCount;
   res.resize(pointsCount);

   for(size_t i = 0; i < pointsCount; i++)
   {
      real px, py, pz, hx, hy, hz;
      
      fin >> _ >> px >> py >> pz >> hx >> hy >> hz;

      res[i] = std::pair<Vector3, Vector3>(Vector3(px, py, pz), Vector3(hx, hy, hz));
   }

   return res;
}

Torus createTorus()
{
   const double torusRadius = 2;
   const double torusSectionWidth = 0.2;
   return Torus(torusRadius, torusSectionWidth, 80, 4, 4);
}

BasisQuadratures readBasisQuadratures()
{
   BasisQuadratures bq;
   std::string bqDir = "basis_quadratures/";

   try
   {
      bq.InitFromTXT(bqDir + "keast 7 x.txt", bqDir + "keast 7 w.txt");
      std::cout << "Quadratures were successfully read. Order = " << bq.order() << std::endl;
   }
   catch(Exeption ex)
   {
      std::cout << ex;
   }

   return bq;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
   for(size_t i = 0; i < vec.size(); i++)
   {
      os << vec[i] << std::endl;
   }

   return os;
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

void comparisonToTelmaIntegrals()
{
   Torus torus = createTorus();
   BasisQuadratures bq = readBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

   auto telmaResults = readTelmaResults("results/telma_results.txt");
   real sumError = 0;
   size_t n = telmaResults.size();

   for(size_t i = 0; i < n; i++)
   {
      auto point = telmaResults[i].first;
      auto telmaB = telmaResults[i].second * math::mu0;

      Vector3 myB = math::calcBioSavartLaplace(current, point, quadratures);
      real error = 100 * (myB - telmaB).length() / telmaB.length();

      std::cout << std::fixed << std::setw(3) << i << " ";
      point.printWithWidth(std::cout, 6);
      std::cout << std::setw(16) << error << std::endl;

      sumError += error;
   }

   std::cout << sumError / n << std::endl;
}

void comparisonToTelmaWithoutTranslation()
{
   Torus torus = createTorus();
   BasisQuadratures bq = readBasisQuadratures();

   auto telmaResults = readTelmaResults("results/telma_results.txt");
   MultipoleSolver multipoleSolver(torus.tetrahedra, bq);
   multipoleSolver.calcLocalMultipolesWithoutTranslation();
   real sumError = 0;
   size_t n = telmaResults.size();

   for(size_t i = 0; i < n; i++)
   {
      auto point = telmaResults[i].first;
      auto telmaB = telmaResults[i].second * math::mu0;

      Vector3 myB = multipoleSolver.calcB(current, point);
      real error = 100 * (myB - telmaB).length() / telmaB.length();

      std::cout << std::fixed << std::setw(3) << i << " ";
      point.printWithWidth(std::cout, 6);
      std::cout << std::setw(16) << error << std::endl;

      sumError += error;
   }

   std::cout << sumError / n << std::endl;
}

void comparisonToTelmaWithTranslation()
{
   Torus torus = createTorus();
   BasisQuadratures bq = readBasisQuadratures();

   auto telmaResults = readTelmaResults("results/telma_results.txt");
   MultipoleSolver multipoleSolver(torus.tetrahedra, bq);
   multipoleSolver.calcLocalMultipolesWithTranslation();
   real sumError = 0;
   size_t n = telmaResults.size();

   for(size_t i = 0; i < n; i++)
   {
      auto point = telmaResults[i].first;
      auto telmaB = telmaResults[i].second * math::mu0;

      Vector3 myB = multipoleSolver.calcB(current, point);
      real error = 100 * (myB - telmaB).length() / telmaB.length();

      std::cout << std::fixed << std::setw(3) << i << " ";
      point.printWithWidth(std::cout, 6);
      std::cout << std::setw(16) << error << std::endl;

      sumError += error;
   }

   std::cout << sumError / n << std::endl;
}

void comparisonBetweenMethodsOnPrecision()
{
   Torus torus = createTorus();
   BasisQuadratures bq = readBasisQuadratures();
   MultipoleSolver multipoleSolver(torus.tetrahedra, bq);
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

   Vector3 point(3, 1, 2);

   Vector3 byIntegration = math::calcBioSavartLaplace(current, point, quadratures);
   
   multipoleSolver.calcLocalMultipolesWithoutTranslation();
   Vector3 byMultipolesWithoutTranslation = multipoleSolver.calcB(current, point);

   multipoleSolver.calcLocalMultipolesWithTranslation();
   Vector3 byMultipolesWithTranslation = multipoleSolver.calcB(current, point);

   std::cout << std::setw(20) << "point ";
   point.printWithWidth(std::cout, 6);
   std::cout << std::scientific << std::endl;
   std::cout << std::setw(40) << "integration " << byIntegration << std::endl;
   std::cout << std::setw(40) << "multipoles w/t translation " << byMultipolesWithoutTranslation << std::endl;
   std::cout << std::setw(40) << "multipoles with translation " << byMultipolesWithTranslation << std::endl;
}

void translationTest()
{
   Vector3 point1(3, 1, 2);
   auto r1 = Harmonics::calcRegularSolidHarmonics(10, point1);
   auto c1 = Harmonics::realToComplex(r1);

   Vector3 point2(4, 5, 1);
   auto r2 = Harmonics::calcRegularSolidHarmonics(10, point1);
   auto c2 = Harmonics::realToComplex(r1);

   auto t =  Harmonics::complexToReal(Harmonics::translate(10, c1, c2));
}

void integrationTimeResearchForMorePoints()
{
   auto torus = createTorus();
   auto bq = readBasisQuadratures();

   for(size_t i = 0; i < 100; i++)
   {

   }
}

void calculationTimeForLocalMultipoles()
{
   auto torus = createTorus();
   auto bq = readBasisQuadratures();

   for(size_t i = 2; i < 15; i++)
   {
      int octreeLeafCapacity = pow(2, i);
      MultipoleSolver multipoleSolver(torus.tetrahedra, bq, octreeLeafCapacity);

      auto start = std::chrono::steady_clock::now();
      multipoleSolver.calcLocalMultipolesWithoutTranslation();
      auto stop = std::chrono::steady_clock::now();
      auto timeWithoutTranslation = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() * 1e-6;

      start = std::chrono::steady_clock::now();
      multipoleSolver.calcLocalMultipolesWithTranslation();
      stop = std::chrono::steady_clock::now();
      auto timeWithTranslation = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() * 1e-6;

      std::cout << octreeLeafCapacity << " " << timeWithoutTranslation << " " << timeWithTranslation << std::endl;
   }
}

int main()
{
   //comparisonToTelmaWithTranslation();
   comparisonBetweenMethodsOnPrecision();
   //calculationTimeForLocalMultipoles();
}
