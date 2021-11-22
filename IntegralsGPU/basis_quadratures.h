#pragma once
#include "real.h"
#include <vector>
#include <string>

namespace triangle_quadratures
{
   class BasisQuadratures
   {
   public:
      std::vector<real> x;
      std::vector<real> y;
      std::vector<real> w;

      int order;

      BasisQuadratures();
      BasisQuadratures(int order);
      void InitFromTXT(std::string coordsFileName, std::string weightsFileName);
   };
};

