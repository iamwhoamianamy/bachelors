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

int main()
{
   const double torusRadius = 2;
   const double torusSectionWidth = 0.2;

   Torus torus(torusRadius, torusSectionWidth, 100, 4, 4);

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
   
   //SphericalHarmonics harmonics(6, Vector3(1, 3, 2));
}
