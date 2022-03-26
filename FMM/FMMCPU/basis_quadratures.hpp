#pragma once
#include "real.hpp"
#include "vector3.hpp"
#include <vector>
#include <string>

class BasisQuadratures
{
private:
   int _order;
   std::vector<Vector3> _values;
   std::vector<real> _w;

public:

   BasisQuadratures();
   void InitFromTXT(std::string coordsFileName, std::string weightsFileame);

   int order() const;
   const std::vector<Vector3>& values() const;
   const std::vector<real>& w() const;

private:
   void InitCoordinatesFromTXT(std::string coordsFileName);
   void InitWeightsFromTXT( std::string weightsFileName);
};