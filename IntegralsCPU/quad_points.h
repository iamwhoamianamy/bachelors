#pragma once
#include <vector>
#include <string>

namespace triangle_quadratures
{
   class QuadPoints
   {
   public:
      std::vector<double> x;
      std::vector<double> y;
      std::vector<double> w;

      int order;

      QuadPoints();
      QuadPoints(int order);
      void InitFromTXT(std::string coordsFileName, std::string weightsFileName);
   };
};

