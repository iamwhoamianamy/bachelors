#include <stdio.h>
#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>

#include "triangle_quadratures.h"
#include "mesh.h"
#include "laplace_solver.h"

using namespace std;
using namespace triangle_quadratures;
using namespace laplace_solver;

const double PI = 3.14159265359;

double f(Vector3 v)
{
   return  1;
}

int main()
{
   Mesh mesh;

   try
   {
      mesh.InitFromOBJ("icoshpere_highres.obj");
   }
   catch(Exeption fileExeption)
   {
      cout << fileExeption;
      exit(1);
   }

   QuadPoints qp(6);
   Vector3 n;

   vector<double> res;
   vector<Vector3> points = { {0.0, 0.0, 0.0} };

   calcIntegralOverMesh(mesh, qp, points, res);

   cout << "True value =       " << u(points[0]) << endl;
   cout << "Calculated value = " << res[0] << endl;

   //cout << "True value =       " << 4.0 * PI * 100 * 100 << endl;
   //cout << "Calculated value = " << calcSurfaceArea(mesh, qp) << endl;

   return 0;
}