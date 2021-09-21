#pragma once
#include <vector>

namespace triangle_quadratures
{
   class QuadPoints
   {
   public:
      std::vector<double> x;
      std::vector<double> y;
      std::vector<double> w;

      int order;

      QuadPoints(int order);
   };
};

