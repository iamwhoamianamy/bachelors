#include "triangle_quadratures.h"
#include <vector>

namespace triangle_quadratures
{
   double calcIntegralOverTriangle(double (*f)(Point), const Triangle& tr, const QuadPoints& qp)
   {
      double result = 0;

      for (size_t i = 0; i < qp.order; i++)
      {
         result += qp.w[i] * f(tr.PointFromST(qp.x[i], qp.y[i]));
      }

      return result * tr.Area();
   }

   void calcIntegralOverMesh(double(*f)(Point),
                             const Mesh& mesh, 
                             const QuadPoints& qp,
                             const vector<Point>& points,
                             std::vector<double>& result)
   {
      
   }
}