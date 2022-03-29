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
   int _n = 10;
   real _eps = 1e-6;
   std::vector<Quadrature> _quadratures;
   Octree* octreeRoot;

public:
   MultipoleSolver(const std::vector<Tetrahedron>& mesh,
                   const BasisQuadratures& basisQuadratures);

   Vector3 calcAFromRoot(real current, const Vector3& point);
   Vector3 calcAWithoutMultipoleTranslation(real current, const Vector3& point);
   Vector3 calcBWithoutMultipoleTranslation(real current, const Vector3& point);
   //Vector3 calcAwithFastMultipoleMethod(real current, const Vector3& point);
   ~MultipoleSolver();

private:
};

