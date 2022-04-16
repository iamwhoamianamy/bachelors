#include <iostream>
#include <chrono>
#include <iomanip>
#include <fstream>

#include "multipole_solver.hpp"
#include "math.hpp"
#include "integration.hpp"
#include "harmonics.hpp"
#include "math.hpp"
#include "translation_algorithms.hpp"

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
   calcLocalMultipolesWithLayersOrMatrices(useGPU, false);
}

void MultipoleSolver::calcLocalMultipolesWithMatrices(bool useGPU)
{
   calcLocalMultipolesWithLayersOrMatrices(useGPU, true);
}

void MultipoleSolver::calcLocalMultipolesWithLayersOrMatrices(bool useGPU, bool useMatrices)
{
   std::vector<std::vector<OctreeNode*>> layers;
   enumerateNodes(octreeRoot, layers, 0);
   calcMultipolesAtLeaves(layers);
   octreeRoot->initAllMultipoleExpansions(n);

   if(_log)
   {
      std::cout << "-------------------------------------------------------" << std::endl;
      std::cout << std::setw(10) << "layer" << std::setw(15) << "mlpl count";
      std::cout << std::setw(15) << "kernel time" << std::setw(15) << "total time" << std::endl;
      std::cout << "-------------------------------------------------------" << std::endl;
      std::cout << std::fixed;
   }

   calcContributionsToHigherLevels(layers, useGPU, useMatrices);
   _multipolesAreReady = true;
}

void MultipoleSolver::enumerateNodes(
   OctreeNode* node,
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

void MultipoleSolver::calcContributionsToHigherLevels(
   const std::vector<std::vector<OctreeNode*>>& layers,
   bool useGPU,
   bool useMatrices)
{
   useMatrices ? calcContributionsToHigherLevelsWithMatrices(layers, useGPU) :
      calcContributionsToHigherLevels(layers, useGPU);
}

void MultipoleSolver::calcMultipolesAtLeaves(
   const std::vector<std::vector<OctreeNode*>>& layers)
{
   for(auto &layer : layers)
   {
      for(auto node : layer)
      {
         if(!node->quadratures().empty())
            node->multipoleExpansion() = math::calcIntegralContribution(
               node->quadratures(), n, node->box().center);
      }
   }
}

void MultipoleSolver::calcContributionsToHigherLevels(
   const std::vector<std::vector<OctreeNode*>>& layers,
   bool useGPU)
{
   for(int l = layers.size() - 1; l >= 1; l--)
   {
      auto start = std::chrono::steady_clock::now();

      if(_log)
      {
         std::cout << std::setw(10) << l << std::setw(15) << layers[l].size();
      }

      std::vector<Vector3> contributions =
         calcContributionsToHigherLevel(layers[l], useGPU);

      for(size_t c = 0; c < layers[l].size(); c++)
      {
         for(size_t j = 0; j < harmonicLength; j++)
         {
            layers[l][c]->parent()->multipoleExpansion().getHarmonic(j) +=
               contributions[c * harmonicLength + j];
         }
      }

      auto stop = std::chrono::steady_clock::now();
      double time = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() * 1e-6;

      if(_log)
      {
         std::cout << std::setw(15) << time << std::endl;
      }
   }
}

std::vector<Vector3> MultipoleSolver::calcContributionsToHigherLevel(
   const std::vector<OctreeNode*>& layer,
   bool useGPU)
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

   auto start = std::chrono::steady_clock::now();

   if(useGPU)
      kernels::translateAllGPU(
         result.data(), regulars.data(), harmonics.data(), layer.size(), n);
   else
      kernels::translateAllCPU(
         result.data(), regulars.data(), harmonics.data(), layer.size(), n);

   auto stop = std::chrono::steady_clock::now();
   double layerTime = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() * 1e-6;

   if(_log)
   {
      std::cout << std::setw(15) << layerTime;
   }

   return result;
}

void MultipoleSolver::calcContributionsToHigherLevelsWithMatrices(
   const std::vector<std::vector<OctreeNode*>>& layers,
   bool useGPU)
{
   for(size_t l = layers.size() - 1; l >= 1; l--)
   {
      auto start = std::chrono::steady_clock::now();

      if(_log)
      {
         std::cout << std::setw(10) << l << std::setw(15) << layers[l].size();
      }

      auto& layer = layers[l];
      auto nodesByOrientation = separateByOrientation(layer);
      auto regularMatrices = calcRegularMatrices(nodesByOrientation);

      double kernelTime = 0;

      for(size_t o = 0; o < 8; o++)
      {
         if(!nodesByOrientation[o].empty())
         {
            auto expansionMatrices =
               getComponentsOfExpansionsInOneOrientation(nodesByOrientation[o]);

            std::vector<ComplexMatrix> translated;
            translated.reserve(3);

            auto kernelStart = std::chrono::steady_clock::now();

            if(useGPU)
            {
               for(size_t c = 0; c < 3; c++)
               {
                  auto regularMatrixAsVector = math::matrixToVector(regularMatrices[o]);
                  auto expansionMatrixAsVector = math::matrixToVector(expansionMatrices[c]);

                  std::vector<std::vector<real>> t(harmonicLength * expansionMatrices[c][0].size());

                  /*kernels::translateAllGPUMatrix(
                     t.data(),
                     regularMatrixAsVector.data(),
                     expansionMatrixAsVector.data(),
                     expansionMatrices[c][0].size(),
                     n);*/

                  //translated.emplace_back(math::vectorToMatrix(t, expansionMatrices[c][0].size()));
               }
            }
            else
            {
               for(size_t c = 0; c < 3; c++)
               {
                  translated.emplace_back(math::mult(
                     regularMatrices[o],
                     expansionMatrices[c]));
               }
            }

            auto kernelStop = std::chrono::steady_clock::now();
            kernelTime += std::chrono::duration_cast<std::chrono::microseconds>
               (kernelStop - kernelStart).count() * 1e-6;

            accountContributions(nodesByOrientation[o], translated);
         }
      }

      auto stop = std::chrono::steady_clock::now();
      double layerTime = std::chrono::duration_cast<std::chrono::microseconds>
         (stop - start).count() * 1e-6;

      if(_log)
      {
         std::cout << std::setw(15) << kernelTime;
         std::cout << std::setw(15) << layerTime << std::endl;
      }
   }
}

Matrix<OctreeNode*> MultipoleSolver::separateByOrientation(
   const std::vector<OctreeNode*>& layer)
{
   Matrix<OctreeNode*> res(8);

   for(auto node : layer)
   {
      for(size_t i = 0; i < 8; i++)
      {
         if(node->parent()->children()[i] == node)
         {
            res[i].push_back(node);
            break;
         }
      }
   }

   return res;
}

std::vector<ComplexMatrix> MultipoleSolver::calcRegularMatrices(
   const Matrix<OctreeNode*>& nodesByOrientation)
{
   std::vector<ComplexMatrix> res;
   res.reserve(8);

   for(size_t i = 0; i < 8; i++)
   {
      auto parent = nodesByOrientation[i][0]->parent();

      if(!nodesByOrientation[i].empty())
      {
         auto translation = nodesByOrientation[i][0]->box().center - parent->box().center;
         auto regularHarmonic = Harmonics::calcRegularSolidHarmonics(n, translation);
         res.emplace_back(matrixFromRegularHarmonic(Harmonics::realToComplex(regularHarmonic)));
      }
   }

   return res;
}

ComplexMatrix MultipoleSolver::matrixFromRegularHarmonic(
   const ComplexHarmonicSeries& regular)
{
   ComplexMatrix res(
      regular.elemCount(),
      std::vector<Complex>(regular.elemCount()));

   for(int l = 0; l <= regular.order(); l++)
   {
      for(int m = -l; m <= l; m++)
      {
         for(int lambda = 0; lambda <= l; lambda++)
         {
            int dl = l - lambda;

            for(int mu = -lambda; mu <= lambda; mu++)
            {
               int dm = m - mu;

               if(-dl <= dm && dm <= dl)
               {
                  res[l * l + l + m][dl * dl + dl + dm] =
                     regular.getHarmonic(lambda * lambda + lambda + mu) *
                     Harmonics::strangeFactor(m, mu);
               }
            }
         }
      }
   }

   return res;
}

std::vector<ComplexMatrix> MultipoleSolver::getComponentsOfExpansionsInOneOrientation(
   const std::vector<OctreeNode*>& nodesByOrientation)
{
   std::vector<ComplexMatrix> res(3, ComplexMatrix(harmonicLength));
   
   for(size_t c = 0; c < 3; c++)
   {
      for(size_t i = 0; i < harmonicLength; i++)
      {
         res[c][i].resize(nodesByOrientation.size());
      }
   }

   for(size_t h = 0; h < nodesByOrientation.size(); h++)
   {
      auto& expansion = nodesByOrientation[h]->multipoleExpansion();

      for(size_t c = 0; c < 3; c++)
      {
         auto ñomplex = Harmonics::realToComplex(Harmonics::separateCoord(expansion, c));

         for(size_t i = 0; i < harmonicLength; i++)
         {
            res[c][i][h] = ñomplex.getHarmonic(i);
         }
      }
   }

   return res;
}

void MultipoleSolver::accountContributions(
   const std::vector<OctreeNode*>& nodesByOrientation,
   const std::vector<ComplexMatrix>& contributions)
{
   for(size_t nodeId = 0; nodeId < nodesByOrientation.size(); nodeId++)
   {
      auto parent = nodesByOrientation[nodeId]->parent();

      auto xComponent = Harmonics::complexToReal(
         ComplexHarmonicSeries(math::getColumn(contributions[0], nodeId)));

      auto yComponent = Harmonics::complexToReal(
         ComplexHarmonicSeries(math::getColumn(contributions[1], nodeId)));

      auto zComponent = Harmonics::complexToReal(
         ComplexHarmonicSeries(math::getColumn(contributions[2], nodeId)));
      
      auto totalContribution = Harmonics::createFormXYZ(
         xComponent, yComponent, zComponent);

      parent->multipoleExpansion().add(totalContribution);
   }
}

void MultipoleSolver::printMatrices(
   const std::vector<ComplexMatrix>& regularMatrices,
   const std::vector<ComplexMatrix>& expansionMatrices)
{
   for(size_t i = 0; i < 8; i++)
   {
      operator<<(std::ofstream("matrices/regular_" + std::to_string(i) +".txt"), 
                 regularMatrices[i]);
   }

   operator<<(std::ofstream("matrices/expansion_x.txt"), expansionMatrices[0]);
   operator<<(std::ofstream("matrices/expansion_y.txt"), expansionMatrices[1]);
   operator<<(std::ofstream("matrices/expansion_z.txt"), expansionMatrices[2]);
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
