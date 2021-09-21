#include "quad_points.h"
#include "exeptions.h"

namespace triangle_quadratures
{
   QuadPoints::QuadPoints(int order) : order(order)
   {
      switch(order)
      {
         case 6:
         {
            x = { 0.659027622374092, 0.659027622374092, 0.231933368553031, 0.231933368553031, 0.109039009072877, 0.109039009072877 };
            y = { 0.231933368553031, 0.109039009072877, 0.659027622374092, 0.109039009072877, 0.659027622374092, 0.231933368553031 };
            w = { 0.166666666666667, 0.166666666666667, 0.166666666666667, 0.166666666666667, 0.166666666666667, 0.166666666666667 };
            break;
         }
         default:
         {
            throw NotImplementedExeption();
         }
      }
   }

};