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
   Mesh sphere;

   try
   {
      sphere.InitFromOBJ("icosphere.obj");
   }
   catch(Exeption fileExeption)
   {
      cout << fileExeption;
      exit(1);
   }

   QuadPoints qp(6);
   Vector3 n;

   vector<double> res;
   vector<Vector3> points = { {0.075, 0.1, 0.025} };

   calcIntegralOverMesh(sphere, qp, points, res);

   cout << "True value =       " << u(points[0]) << endl;
   cout << "Calculated value = " << res[0] << endl;

   /*cout << "True value =       " << 4.0 * PI << endl;
   cout << "Calculated value = " << calcIntegralOverMesh(f, sphere, qp) << endl;*/

   

   return 0;
}