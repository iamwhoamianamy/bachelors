#pragma once
#include "real.hpp"
#include <vector>
#include <string>

class BasisQuadratures
{
public:
   std::vector<real> x;
   std::vector<real> y;
   std::vector<real> z;
   std::vector<real> w;

   int order;

   BasisQuadratures();
   void InitFromTXT(std::string coordsFileName, std::string weightsFileName);
   void InitCoordinatesFromTXT(std::string coordsFileName);
   void InitWeightsFromTXT( std::string weightsFileName);
};