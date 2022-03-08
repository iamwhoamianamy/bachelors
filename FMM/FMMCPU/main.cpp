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

using namespace std;
const double PI = 3.14159265359;

int main()
{
   /*const double torusRadius = 3;
   const double torusSectionWidth = 1;

   Torus torus(torusRadius, torusSectionWidth, 20, 2, 4);

   BasisQuadratures bq;
   string bqDir = "basis_quadratures/";

   try
   {
      bq.InitFromTXT(bqDir + "keast 7 x.txt", bqDir + "keast 7 w.txt");
      cout << "Quadratures were successfully read. Order = " << bq.order << endl;
   }
   catch(Exeption ex)
   {
      cout << ex;
   }*/

   //Vector3 a(0, 0, 0);
   //Vector3 b(0, 1, 0);
   //Vector3 c(1, 0, 0);
   //Vector3 d(0, 0, 1);

   //Tetrahedron tetr(a, b, c, d);

   //cout << tetr.Volume();

   /*double res = 0;

   for(auto &t :torus.tetrahedra)
   {
      res += t.Volume();
   }

   double majorTorusRadius = PI * pow(torusRadius + torusSectionWidth / 2, 2);
   double minorTorusRadius = PI * pow(torusRadius - torusSectionWidth / 2, 2);

   double trueTorusVolume = (majorTorusRadius - minorTorusRadius) * torusSectionWidth;

   cout << "True torus volume = " << trueTorusVolume << endl;
   cout << "Generated torus volume = " << res << endl;*/

   //auto vec = Vector3(200, 200, 0);
   //std::vector<Vector3> vectors;
   //vectors.push_back(vec);

   std::cout << math::calcFactorial(3);
}
