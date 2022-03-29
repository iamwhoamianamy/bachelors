#include <iostream>
#include <string>
#include <fstream>

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

void runCalculations()
{
   const double torusRadius = 2;
   const double torusSectionWidth = 0.2;

   Torus torus(torusRadius, torusSectionWidth, 20, 4, 4);

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

   int n = 10;
   real current = 5;
   Vector3 point(2, 1, 3);

   MultipoleSolver multipoleSolver(torus.tetrahedra, bq);

   Vector3 aIntegrRes = calcAViaSimpleIntegration(current, point, torus.tetrahedra, bq);
   Vector3 aMultRes = multipoleSolver.calcAWithoutMultipoleTranslation(current, point);

   cout << fixed;
   cout << "A by simple integration:        " << aIntegrRes << endl;
   cout << "A by multipole method:          " << aMultRes << endl;

   Vector3 hIntegrRes = calcBioSavartLaplace(current, point, torus.tetrahedra, bq) / math::mu0;
   Vector3 hMultRes = multipoleSolver.calcBWithoutMultipoleTranslation(current, point) / math::mu0;
   cout << endl;
   cout << "H by simple integration:        " << hIntegrRes << endl;
   cout << "H by multipole method:          " << hMultRes << endl;
}

void runTest()
{
   real current = 5;
   Vector3 point(1, 2, 3);
   auto solidHarmonics = Harmonics::calcSolidHarmonics(15, point, false);
}

int main()
{   
   runCalculations();
}
