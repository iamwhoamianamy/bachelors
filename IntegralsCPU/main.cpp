#include <stdio.h>
#include <iostream>
#include <cmath>
#include <iomanip>

#include "triangle_quadratures.h"
#include "mesh.h"

using namespace std;
using namespace triangle_quadratures;

const double PI = 3.14159265359;

double f(Point p)
{
   return  1;
}

int main()
{
   Mesh sphere("icosphere.txt");
   QuadPoints qp = QuadPoints(6);

   double sum = 0;
   double answer = 12.5663706144;

   cout << "i:" << setw(4) << " area:" << setw(14) << "integral" << setw(14) << endl;

   for(size_t i = 0; i < sphere.TriangleCount(); i++)
   {
      Triangle tr = sphere.GetTriangle(i);
      double area = tr.Area();
      double integral = calcIntegral(f, tr, qp);
      cout << setw(4) << i << setw(14) << area << setw(14) << integral << endl;
      sum += integral;
   }

   cout << "Sum is " << sum << endl;

   return 0;
}