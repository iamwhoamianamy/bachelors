#include <stdio.h>
#include <iostream>

#include "qwer.h"

using namespace std;

double f(double x, double y)
{
   return y*y;
}

int main()
{
   QuadPoints qp = QuadPoints(6);
   Triangle tr = Triangle(0, 1,
                          2, -1,
                          2, 3);

   cout << "Area is " << tr.GetArea() << endl;
   cout << "Integral is " << integral(f, tr, qp) << endl;

   return 0;
}