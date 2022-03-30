#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

#include "vector3.hpp"
#include "tetrahedron.hpp"
#include "hexahedron.hpp"
#include "torus.hpp"
#include "basis_quadratures.hpp"
#include "exeptions.hpp"
#include "math.hpp"
#include "spherical_harmonics.hpp"
#include "multipole_solver.hpp"

using namespace math;

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
   return Torus(torusRadius, torusSectionWidth, 40, 4, 4);
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

void comparisonToTelmaWithoutTranslation()
{
   Torus torus = createTorus();
   BasisQuadratures bq = readBasisQuadratures();

   auto telmaResults = readTelmaResults("results/telma_results.txt");
   real current = 5;
   MultipoleSolver multipoleSolver(torus.tetrahedra, bq);
   multipoleSolver.calcLocalMultipolesWithoutTranslation();
   real sumErrorForIntegrals = 0;
   real sumErrorForMultipoles = 0;
   size_t n = telmaResults.size();

   for(size_t i = 0; i < n; i++)
   {
      auto point = telmaResults[i].first;
      auto telmaB = telmaResults[i].second * math::mu0;

      Vector3 myBIntegrals = math::calcBioSavartLaplace(current, point, torus.tetrahedra, bq);
      Vector3 myBMultipoles = multipoleSolver.calcB(current, point);
      real errorForIntegrals = 100 * (myBIntegrals - telmaB).length() / telmaB.length();
      real errorForMultipoles = 100 * (myBMultipoles - telmaB).length() / telmaB.length();

      std::cout << std::fixed << std::setw(3) << i << " ";
      point.printWithWidth(std::cout, 6);
      std::cout << std::scientific << std::setw(16) << errorForIntegrals << " ";
      std::cout << std::setw(16) << errorForMultipoles << std::endl;

      sumErrorForIntegrals += errorForIntegrals;
      sumErrorForMultipoles += errorForMultipoles;
   }

   std::cout << sumErrorForIntegrals / n << " " << sumErrorForMultipoles / n << std::endl;
}

void comparisonBetweenMethods()
{
   Torus torus = createTorus();
   BasisQuadratures bq = readBasisQuadratures();
   MultipoleSolver multipoleSolver(torus.tetrahedra, bq);
   real current = 5;
   Vector3 point(3, 1, 2);

   Vector3 byIntegration = math::calcBioSavartLaplace(current, point, torus.tetrahedra, bq);
   
   multipoleSolver.calcLocalMultipolesWithoutTranslation();
   Vector3 byMultipolesWithoutTranslation = multipoleSolver.calcB(current, point);

   multipoleSolver.calcLocalMultipolesWithTranslation();
   Vector3 byMultipolesWithTranslation = multipoleSolver.calcB(current, point);

   std::cout << std::setw(20) << "point ";
   point.printWithWidth(std::cout, 6);
   std::cout << std::scientific << std::endl;
   std::cout << std::setw(40) << "integration " << byIntegration << std::endl;
   std::cout << std::setw(40) << "multipoles w/t translation " << byMultipolesWithoutTranslation << std::endl;
   std::cout << std::setw(40) << "multipoles with translation " << byMultipolesWithTranslation;
}

void test()
{
   Vector3 point1(3, 1, 2);
   auto r1 = Harmonics::calcRegularSolidHarmonics(10, point1);
   auto c1 = Harmonics::realToComplex(r1);

   Vector3 point2(4, 5, 1);
   auto r2 = Harmonics::calcRegularSolidHarmonics(10, point1);
   auto c2 = Harmonics::realToComplex(r1);

   auto t =  Harmonics::complexToReal(Harmonics::translate(10, c1, c2));
}

int main()
{
   comparisonBetweenMethods();
}
