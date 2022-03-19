#pragma once
#include "real.hpp"
#include <vector>
#include <string>

class BasisQuadratures
{
private:
   int _order;
   std::vector<real> _x;
   std::vector<real> _y;
   std::vector<real> _z;
   std::vector<real> _w;

public:

   BasisQuadratures();
   void InitFromTXT(std::string coordsFileName, std::string weightsFileame);

   int order() const;
   const std::vector<real>& x() const;
   const std::vector<real>& y() const;
   const std::vector<real>& z() const;
   const std::vector<real>& w() const;

private:
   void InitCoordinatesFromTXT(std::string coordsFileName);
   void InitWeightsFromTXT( std::string weightsFileName);
};