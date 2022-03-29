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

using namespace std;
using namespace math;

vector<pair<Vector3, Vector3>> readTelmaResults(const string& filename)
{
   vector<pair<Vector3, Vector3>> res;
   ifstream fin(filename);
   size_t pointsCount;
   string _;
   fin >> _ >> _ >> _;
   fin >> pointsCount;
   res.resize(pointsCount);

   for(size_t i = 0; i < pointsCount; i++)
   {
      real px, py, pz, hx, hy, hz;
      
      fin >> _ >> px >> py >> pz >> hx >> hy >> hz;

      res[i] = pair<Vector3, Vector3>(Vector3(px, py, pz), Vector3(hx, hy, hz));
   }

   return res;
}

void runCalculations()
{
   const double torusRadius = 2;
   const double torusSectionWidth = 0.2;

   Torus torus(torusRadius, torusSectionWidth, 40, 4, 4);

   BasisQuadratures bq;
   string bqDir = "basis_quadratures/";

   try
   {
      bq.InitFromTXT(bqDir + "keast 7 x.txt", bqDir + "keast 7 w.txt");
      cout << "Quadratures were successfully read. Order = " << bq.order() << endl;
   }
   catch(Exeption ex)
   {
      cout << ex;
   }

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
      Vector3 myBMultipoles = multipoleSolver.calcBWithoutMultipoleTranslation(current, point);
      real errorForIntegrals = 100 * (myBIntegrals - telmaB).length() / telmaB.length();
      real errorForMultipoles = 100 * (myBMultipoles - telmaB).length() / telmaB.length();

      cout << fixed << setw(3) << i << " ";
      point.printWithWidth(cout, 6);
      cout << scientific << setw(16) << errorForIntegrals << " ";
      cout << setw(16) << errorForMultipoles << endl;

      sumErrorForIntegrals += errorForIntegrals;
      sumErrorForMultipoles += errorForMultipoles;
   }

   cout << sumErrorForIntegrals / n << " " << sumErrorForMultipoles / n << endl;
}

int main()
{
   runCalculations();
}
