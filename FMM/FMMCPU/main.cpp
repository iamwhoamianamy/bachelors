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

using namespace std;
using namespace math;

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

   real current = 5;
   Vector3 point(2, 1, 3);
   Vector3 H = calcBioSavartLaplace(current, point, torus.tetrahedra, bq);

   /*for(size_t i = 0; i < torus.tetrahedra.size(); i++)
   {
      cout << i << " " << torus.tetrahedra[i].volume() << endl;
   }*/

   cout << endl;
   cout << H / 1.256e-6 << endl;
}

void runTest()
{
   Vector3 point(1, 2, 3);
   auto solidHarmonics = SphericalHarmonics::calcSolidHarmonics(15, point, true);
}

int main()
{   
   runTest();
}
