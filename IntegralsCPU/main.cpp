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
   return v.x + v.y + v.z;
}

int main()
{
   Mesh mesh;

   try
   {
      mesh.InitFromOBJ("meshes/icosphere.obj");
      //mesh.InitFromOBJ("meshes/cube_highres.obj");
   }
   catch(Exeption fileExeption)
   {
      cout << fileExeption;
      exit(1);
   }

   QuadPoints qp;
  
   try
   {
      qp.InitFromTXT("quadratures/gauss15_xy.txt", "quadratures/gauss15_w.txt");
   }
   catch(Exeption fileExeption)
   {
      cout << fileExeption;
      exit(1);
   }

   Vector3 n;

   vector<double> res;
   vector<Vector3> points = { {0.8, 0.20, 0.00} };

   calcIntegralOverMesh(mesh, qp, points, res);

   double true_value = u(points[0]);
   double calc_value = res[0];
   double error = abs((true_value - calc_value) / true_value);

   cout << "Integral:" << endl;
   cout << fixed;
   cout << "True value =       " << setw(16) << true_value << endl;
   cout << "Calculated value = " << setw(16) << calc_value << endl;
   cout << scientific;
   cout << "Error            = " << setw(16) << error << endl;

   //true_value = 4.0 * PI;
   ////true_value = 4.0 * 6;
   //calc_value = calcSurfaceArea(mesh, qp);
   //error = abs(true_value - calc_value);

   //cout << endl << "Surface:" << endl;
   //cout << fixed;
   //cout << "True value =       " << setw(16) << true_value << endl;
   //cout << "Calculated value = " << setw(16) << calc_value << endl;
   //cout << scientific;
   //cout << "Error            = " << setw(16) << error << endl;

   //cout << calcIntegralOverMesh(f, mesh, qp) << endl;
   //cout << calcSurfaceArea(mesh, qp);

   return 0;
}