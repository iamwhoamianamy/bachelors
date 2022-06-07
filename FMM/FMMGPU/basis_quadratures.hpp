#pragma once
#include "real.hpp"
#include "vector3.cuh"
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
   void initFromTXT(std::string coordsFileName, std::string weightsFileName);

   int order() const;
   const Vector3& values(size_t i) const;
   real w(size_t i) const;

private:
   void initCoordinatesFromTXT(std::string coordsFileName);
   void initWeightsFromTXT( std::string weightsFileName);
};