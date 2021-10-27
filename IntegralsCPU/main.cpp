#include <stdio.h>
#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>
#include <chrono>

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
      mesh.InitFromOBJ("../meshes/icosphere_highres.obj");
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
      qp.InitFromTXT("../quadratures/gauss15_xy.txt", "../quadratures/gauss15_w.txt");
   }
   catch(Exeption fileExeption)
   {
      cout << fileExeption;
      exit(1);
   }

   Vector3 n;

   const int points_count = 1000;
   vector<double> res;
   vector<Vector3> points(points_count);

   for(size_t i = 0; i < points_count; i++)
   {
      points[i] = { 0.8 / points_count * (i + 1), 0.20, 0.00 };
   }

   auto start = std::chrono::steady_clock::now();

   calcIntegralOverMesh(mesh, qp, points, res);

   auto stop = std::chrono::steady_clock::now();

   cout << "Calculation time: " << chrono::duration_cast<chrono::microseconds>(stop - start).count() * 1e-6 << endl << endl;

   for(size_t i = 0; i < points_count; i++)
   {
      cout << "Point: " << scientific << points[i].x << " " << points[i].y << " " << points[i].z << endl;

      double true_value = u(points[i]);
      double calc_value = res[i];
      double error = abs((true_value - calc_value) / true_value);

      cout << "Integral:" << endl;
      cout << fixed;
      cout << "True value =       " << setw(16) << true_value << endl;
      cout << "Calculated value = " << setw(16) << calc_value << endl;
      cout << scientific;
      cout << "Error            = " << setw(16) << error << endl;
   }
   

   return 0;
}