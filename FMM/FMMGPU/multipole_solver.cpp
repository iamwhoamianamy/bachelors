#include <iostream>
#include <chrono>
#include <iomanip>

#include "multipole_solver.hpp"
#include "math.hpp"
#include "integration.hpp"
#include "harmonics.hpp"
#include "math.hpp"
#include "kernel_callers.hpp"

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

void MultipoleSolver::calcLocalMultipolesWithLayers(bool useGPU)
{
   std::vector<std::vector<OctreeNode*>> layers;
   enumerateNodes(octreeRoot, layers, 0);
   calcMultipolesAtLeaves(layers);
   octreeRoot->initAllMultipoleExpansions(n);
   calcContributionsToHigherLevel(layers, useGPU);
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

void MultipoleSolver::calcContributionsToHigherLevel(
   const std::vector<std::vector<OctreeNode*>>& layers, bool useGPU)
{
   std::cout << "-----------------------------------------------" << std::endl;
   std::cout << std::setw(20) << "layer" << std::setw(15) << "mlpl count";
   std::cout << std::setw(10) << "time" << std::endl;
   std::cout << "-----------------------------------------------" << std::endl;
   std::cout << std::fixed;

   for(int i = layers.size() - 1; i >= 1; i--)
   {
      auto start = std::chrono::steady_clock::now();
      std::vector<Vector3> contributions =
         calcContributionsToHigherLevel(layers[i], useGPU);
      auto stop = std::chrono::steady_clock::now();
      auto time = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() * 1e-6;

      std::cout << std::setw(20) << i << std::setw(15) << layers[i].size();
      std::cout << std::setw(10) << time << std::endl;

      for(size_t c = 0; c < layers[i].size(); c++)
      {
         for(size_t j = 0; j < harmonicLength; j++)
         {
            layers[i][c]->parent()->multipoleExpansion().getHarmonic(j) +=
               contributions[c * harmonicLength + j];
         }
      }
   }
}

std::vector<Vector3> MultipoleSolver::calcContributionsToHigherLevel(
   const std::vector<OctreeNode*>& layer, bool useGPU)
{
   std::vector<Vector3> harmonics(layer.size() * harmonicLength);
   std::vector<real> regulars(layer.size() * harmonicLength);

   for(size_t i = 0; i < layer.size(); i++)
   {
      Vector3 translation = layer[i]->box().center - layer[i]->parent()->box().center;
      auto regular = Harmonics::calcRegularSolidHarmonics(n, translation);

      for(size_t j = 0; j < harmonicLength; j++)
      {
         regulars[i * harmonicLength + j] = regular.getHarmonic(j);
         harmonics[i * harmonicLength + j] = layer[i]->multipoleExpansion().getHarmonic(j);
      }
   }

   std::vector<Vector3> result(layer.size() * harmonicLength);

   if(useGPU)
      kernels::translateAllGPU(result.data(), regulars.data(), harmonics.data(), layer.size(), n);
   else
      kernels::translateAllCPU(result.data(), regulars.data(), harmonics.data(), layer.size(), n);

   return result;
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
