#include <iostream>
#include "multipole_solver.hpp"
#include "math.hpp"
#include "integration.hpp"
#include "harmonics.hpp"
#include "math.hpp"

MultipoleSolver::MultipoleSolver(std::vector<Quadrature>& quadratures,
                                 size_t octreeLeafCapacity) :
   _quadratures(quadratures), octreeLeafCapacity(octreeLeafCapacity)
{
   octreeRoot = new OctreeNode(Box(Vector3(0, 0, 0), Vector3(3, 3, 3)), octreeLeafCapacity);
   octreeRoot->insert(_quadratures);
}

void MultipoleSolver::calcLocalMultipolesWithoutTranslation()
{
   octreeRoot->calcLocalMultipolesWithoutTranslation(n);
   _multipolesAreReady = true;
}

void MultipoleSolver::calcLocalMultipolesWithComplexTranslation()
{
   octreeRoot->calcLocalMultipolesWithComplexTranslation(n);
   _multipolesAreReady = true;
}

void MultipoleSolver::calcLocalMultipolesWithRealTranslation()
{
   octreeRoot->calcLocalMultipolesWithRealTranslation(n);
   _multipolesAreReady = true;
}

void MultipoleSolver::calcLocalMultipolesWithLayers()
{
   std::vector<std::vector<OctreeNode*>> layers;
   enumerateNodes(octreeRoot, layers, 0);
   calcMultipolesAtLeaves(layers);
   octreeRoot->initAllMultipoleExpansions(n);

   for(int i = layers.size() - 1; i >= 1; i--)
   {
      auto contributions = calcContributionsToHigherLevel(layers[i]);

      for(size_t c = 0; c < contributions.size(); c++)
      {
         layers[i][c]->parent()->multipoleExpansion().add(contributions[c]);
      }
   }

   _multipolesAreReady = true;
}

void MultipoleSolver::enumerateNodes(OctreeNode* node,
                                     std::vector<std::vector<OctreeNode*>>& layers,
                                     size_t currentLayerId)
{
   if(!node->quadratures().empty() || node->isSubdivided())
   {
      if(layers.size() <= currentLayerId)
         layers.push_back(std::vector<OctreeNode*>());

      layers[currentLayerId].push_back(node);

      for(auto child : node->children())
      {
         enumerateNodes(child, layers, currentLayerId + 1);
      }
   }
}

void MultipoleSolver::calcMultipolesAtLeaves(
   const  std::vector<std::vector<OctreeNode*>>& layers)
{
   for(auto layer : layers)
   {
      for(auto node : layer)
      {
         if(!node->quadratures().empty())
            node->multipoleExpansion() = math::calcIntegralContribution(
               node->quadratures(), n, node->box().center);
      }
   }
}

std::vector<HarmonicSeries<Vector3>> MultipoleSolver::calcContributionsToHigherLevel(
   const std::vector<OctreeNode*>& layer)
{
   std::vector<HarmonicSeries<Vector3>> res;

   for(auto node : layer)
   {
      res.push_back(Harmonics::translateWithReal(
         node->multipoleExpansion(),
         node->box().center - node->parent()->box().center));
   }

   return res;
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
