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
   const Vector3& values(size_t i) const;
   const real w(size_t i) const;

private:
   void InitCoordinatesFromTXT(std::string coordsFileName);
   void InitWeightsFromTXT( std::string weightsFileName);
};