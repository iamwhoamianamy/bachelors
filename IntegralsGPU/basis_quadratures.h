#pragma once
#include <vector>
#include <string>

namespace triangle_quadratures
{
   class BasisQuadratures
   {
   public:
      std::vector<float> x;
      std::vector<float> y;
      std::vector<float> w;

      int order;

      BasisQuadratures();
      BasisQuadratures(int order);
      void InitFromTXT(std::string coordsFileName, std::string weightsFileName);
   };
};

