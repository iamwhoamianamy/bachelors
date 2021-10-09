#include <stdio.h>
#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>

#include "triangle_quadratures.h"
#include "mesh.h"

using namespace std;
using namespace triangle_quadratures;

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

   calcIntegralOverMesh(f, sphere, qp, { 0, 0, 0 }, res);

   cout << res[0];

   return 0;
}