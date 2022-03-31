#pragma once
#include <vector>
#include "vector3.hpp"
#include "tetrahedron.hpp"
#include "basis_quadratures.hpp"
#include "quadrature.hpp"
#include "octree.hpp"

class MultipoleSolver
{
private:
   std::vector<Quadrature>& _quadratures;
   Octree* octreeRoot;
   bool _multipolesAreReady = false;

public:
   const int n = 10;
   const real eps = 1e-6;
   const size_t octreeLeafCapacity;

   MultipoleSolver(std::vector<Quadrature>& quadratures,
                   size_t octreeLeafCapacity = 1000);

   void calcLocalMultipolesWithoutTranslation();
   void calcLocalMultipolesWithTranslation();
   Vector3 calcA(real current, const Vector3& point);
   Vector3 calcB(real current, const Vector3& point);
   ~MultipoleSolver();

private:
};

