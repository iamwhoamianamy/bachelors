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
   Vector3 point(2, 1.6, 0);


   /*Vector3 H = calcBioSavartLaplace(current, point, torus.tetrahedra, bq);
   cout << endl;
   cout << H / 1.256e-6 << endl;*/

   MultipoleSolver multipoleSolver(torus.tetrahedra, bq);

   Vector3 integrRes = calcAViaSimpleIntegration(current, point, torus.tetrahedra, bq);
   Vector3 multRes = calcAViaMultipoleMethod(current, point, torus.tetrahedra, bq, n);
   Vector3 treeRes = multipoleSolver.calcAWithMultipoleMethod(current, point);

   cout << fixed;
   cout << "Simple integration:        " << integrRes << endl;
   cout << "Multipole method:          " << multRes << endl;
   cout << "Multipole with tree:       " << treeRes << endl;
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
