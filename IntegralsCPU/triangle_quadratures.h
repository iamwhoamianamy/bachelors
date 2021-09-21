#pragma once
#include <vector>
#include "triangle.h"
#include "exeptions.h"
#include "quad_points.h"
#include "point.h"


using namespace std;


namespace triangle_quadratures
{
   double calcIntegral(double (*f)(Point), const Triangle& tr, const QuadPoints& qp)
   {
      const int order = qp.order;
      double result = 0;

      for(size_t i = 0; i < order; i++)
      {
         result += qp.w[i] * f(tr.PointFromST(qp.x[i], qp.y[i]));
      }

      return result * tr.Area();
   }
}