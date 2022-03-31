#include <iostream>
#include "multipole_solver.hpp"
#include "math.hpp"
#include "harmonics.hpp"
#include "math.hpp"

MultipoleSolver::MultipoleSolver(std::vector<Quadrature>& quadratures,
                                 size_t octreeLeafCapacity) :
   _quadratures(quadratures), octreeLeafCapacity(octreeLeafCapacity)
{
   octreeRoot = new Octree(Box(Vector3(0, 0, 0), Vector3(3, 3, 3)), octreeLeafCapacity);
   octreeRoot->insert(_quadratures);
}

void MultipoleSolver::calcLocalMultipolesWithoutTranslation()
{
   octreeRoot->calcLocalMultipolesWithoutTranslation(n);
   _multipolesAreReady = true;
}

void MultipoleSolver::calcLocalMultipolesWithTranslation()
{
   octreeRoot->calcLocalMultipolesWithTranslation(n);
   _multipolesAreReady = true;
}

Vector3 MultipoleSolver::calcA(real current, const Vector3& point)
{
   if(!_multipolesAreReady)
      throw new std::exception("Multipoles are not ready!");

   return octreeRoot->calcA(point) / (4.0 * math::PI) * current;
}

Vector3 MultipoleSolver::calcB(real current, const Vector3& point)
{
   if(!_multipolesAreReady)
      throw new std::exception("Multipoles are not ready!");

   return octreeRoot->caclRot(point) / (4.0 * math::PI) * current * math::mu0;
}

MultipoleSolver::~MultipoleSolver()
{
   delete octreeRoot;
}
