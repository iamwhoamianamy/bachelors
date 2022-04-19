#include <iostream>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <omp.h>

#include "multipole_solver.hpp"
#include "math.hpp"
#include "integration.hpp"
#include "harmonics.hpp"
#include "math.hpp"
#include "translation_algorithms.hpp"
#include "kernels.cuh"
#include "testing_helpers.hpp"

MultipoleSolver::MultipoleSolver(std::vector<Quadrature>& quadratures,
                                 size_t octreeLeafCapacity) :
   _quadratures(quadratures), octreeLeafCapacity(octreeLeafCapacity)
{
   octreeRoot = new OctreeNode(Box(Vector3(0, 0, 0), Vector3(3, 3, 3)), octreeLeafCapacity);
   octreeRoot->insert(_quadratures);
}

void MultipoleSolver::calcLocalMultipolesWithoutTranslation()
{
   octreeRoot->calcLocalMultipolesWithoutTranslation(harmonicOrder);
   _multipolesAreReady = true;
}

void MultipoleSolver::calcLocalMultipolesWithComplexTranslation()
{
   octreeRoot->calcLocalMultipolesWithComplexTranslation(harmonicOrder);
   _multipolesAreReady = true;
}

void MultipoleSolver::calcLocalMultipolesWithRealTranslation()
{
   octreeRoot->calcLocalMultipolesWithRealTranslation(harmonicOrder);
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

void MultipoleSolver::calcLocalMultipolesWithLayersOrMatrices(
   bool useGPU,
   bool useMatrices)
{
   std::vector<std::vector<OctreeNode*>> layers;
   enumerateNodes(octreeRoot, layers, 0);
   calcMultipolesAtLeaves(layers);
   octreeRoot->initAllMultipoleExpansions(harmonicOrder);

   if(_log)
   {
      if(useMatrices)
         std::cout << "-----------------------matrices------------------------" << std::endl;
      else
         std::cout << "------------------------layers-------------------------" << std::endl;

      std::cout << std::setw(10) << "layer" << std::setw(15) << "mlpl count";
      std::cout << std::setw(15) << "kernel time" << std::setw(15) << "total time" << std::endl;
      std::cout << "-------------------------------------------------------" << std::endl;
      std::cout << std::fixed;
   }

   calcContributionsToHigherLayers(layers, useGPU, useMatrices);
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

void MultipoleSolver::calcContributionsToHigherLayers(
   const std::vector<std::vector<OctreeNode*>>& layers,
   bool useGPU,
   bool useMatrices)
{
   useMatrices ? calcContributionsToHigherLevelsWithMatrices(layers, useGPU) :
      calcContributionsToHigherLayers(layers, useGPU);
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
               node->quadratures(), harmonicOrder, node->box().center);
      }
   }
}

void MultipoleSolver::calcContributionsToHigherLayers(
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
         calcContributionsToHigherLayer(layers[l], useGPU);

      for(size_t c = 0; c < layers[l].size(); c++)
      {
         for(int h = 0; h < harmonicLength; h++)
         {
            layers[l][c]->parent()->multipoleExpansion().getHarmonic(h) +=
               contributions[c * harmonicLength + h];
         }
      }

      auto stop = std::chrono::steady_clock::now();
      double time = test::getTime(start, stop);

      if(_log)
      {
         std::cout << std::setw(15) << time << std::endl;
      }
   }
}

std::vector<Vector3> MultipoleSolver::calcContributionsToHigherLayer(
   const std::vector<OctreeNode*>& layer,
   bool useGPU)
{
   std::vector<Vector3> harmonics(layer.size() * harmonicLength);
   std::vector<real> regulars(layer.size() * harmonicLength);

   for(size_t i = 0; i < layer.size(); i++)
   {
      Vector3 translation = layer[i]->box().center - layer[i]->parent()->box().center;
      auto regular = Harmonics::calcRegularSolidHarmonics(harmonicOrder, translation);

      std::copy(regular.data().begin(), 
                regular.data().end(),
                regulars.begin() + i * harmonicLength);

      std::copy(layer[i]->multipoleExpansion().data().begin(), 
                layer[i]->multipoleExpansion().data().end(),
                harmonics.begin() + i * harmonicLength);
   }

   std::vector<Vector3> result(layer.size() * harmonicLength);

   auto start = std::chrono::steady_clock::now();

   if(useGPU)
   {
      kernels::translateAllGPU(
         result.data(),
         regulars.data(),
         harmonics.data(),
         layer.size(),
         harmonicOrder);
   }
   else
   {
      kernels::translateAllCPU(
         result.data(),
         regulars.data(),
         harmonics.data(),
         layer.size(),
         harmonicOrder);
   }

   auto stop = std::chrono::steady_clock::now();
   double layerTime = test::getTime(start, stop);

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
      double kernelTime = 0;

      auto start = std::chrono::steady_clock::now();

      if(_log)
      {
         std::cout << std::setw(10) << l << std::setw(15) << layers[l].size();
      }

      auto& layer = layers[l];
      auto nodesByOrientation = separateNodesByOrientation(layer);
      auto regularVectors = calcRegularMatricesForLayerAsVectors(
         nodesByOrientation);
      
      for(size_t o = 0; o < 8; o++)
      {
         if(!nodesByOrientation[o].empty())
         {
            size_t nodesCount = nodesByOrientation[o].size();

            auto expansionVectors =
               getExpansionsInOneOrientationAsVectors(
                  nodesByOrientation[o]);

            ComplexMatrix translated;
            translated.reserve(3);

            for(size_t c = 0; c < 3; c++)
            {
               std::vector<Complex> t(harmonicLength * nodesCount);

               auto kernelStart = std::chrono::steady_clock::now();

               if(useGPU)
               {
                  kernels::translateAllGPUMatrixCuBLAS(
                     t.data(),
                     expansionVectors[c].data(),
                     regularVectors[o].data(),
                     nodesCount,
                     harmonicOrder);
               }
               else
               {
                  kernels::translateAllCPUMatrix(
                     t,
                     expansionVectors[c],
                     regularVectors[o],
                     nodesCount,
                     harmonicOrder);
                  /*kernels::translateAllCPUMatrixBLAS(
                     t.data(),
                     expansionVectors[c].data(),
                     regularVectors[o].data(),
                     nodesCount,
                     harmonicOrder);*/
               }

               auto kernelStop = std::chrono::steady_clock::now();
               kernelTime += std::chrono::duration_cast<std::chrono::microseconds>
                  (kernelStop - kernelStart).count() * 1e-6;

               translated.emplace_back(t);
            }

            accountChildrenContributions(
               nodesByOrientation[o],
               translated);
         }
      }

      auto stop = std::chrono::steady_clock::now();
      double layerTime = test::getTime(start, stop);

      if(_log)
      {
         std::cout << std::setw(15) << kernelTime;
         std::cout << std::setw(15) << layerTime << std::endl;
      }
   }
}

Matrix<OctreeNode*> MultipoleSolver::separateNodesByOrientation(
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

std::vector<ComplexMatrix> MultipoleSolver::calcRegularMatricesForLayer(
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

         auto regularHarmonics = Harmonics::calcRegularSolidHarmonics(
            harmonicOrder,
            translation);

         res.emplace_back(formMatrixFromRegularHarmonics(
            Harmonics::realToComplex(regularHarmonics)));
      }
   }

   return res;
}

ComplexMatrix MultipoleSolver::calcRegularMatricesForLayerAsVectors(
   const Matrix<OctreeNode*>& nodesByOrientation,
   size_t padding)
{
   ComplexMatrix res;
   res.reserve(8);

   for(int i = 0; i < 8; i++)
   {
      auto parent = nodesByOrientation[i][0]->parent();

      if(!nodesByOrientation[i].empty())
      {
         auto translation = nodesByOrientation[i][0]->box().center - parent->box().center;

         auto regularHarmonics = Harmonics::calcRegularSolidHarmonics(
            harmonicOrder, translation);

         res.emplace_back(formMatrixFromRegularHarmonicsAsVectors(
            Harmonics::realToComplex(regularHarmonics), padding));
      }
   }

   return res;
}

ComplexMatrix MultipoleSolver::formMatrixFromRegularHarmonics(
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
                  res[dl * dl + dl + dm][l * l + l + m] =
                     regular.getHarmonic(lambda * lambda + lambda + mu) *
                     Harmonics::strangeFactor(m, mu);
               }
            }
         }
      }
   }

   return res;
}

std::vector<Complex> MultipoleSolver::formMatrixFromRegularHarmonicsAsVectors(
   const ComplexHarmonicSeries& regular,
   size_t padding)
{
   size_t width = math::nextDevisible(
      harmonicLength,
      padding);

   std::vector<Complex> res(width * width);

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
                  res[(l * l + l + m) + (dl * dl + dl + dm) * width] =
                     regular.getHarmonic(lambda * lambda + lambda + mu) *
                     Harmonics::strangeFactor(m, mu);
               }
            }
         }
      }
   }

   return res;
}

std::vector<ComplexMatrix> MultipoleSolver::getExpansionsInOneOrientation(
   const std::vector<OctreeNode*>& nodesByOrientation)
{
   size_t nodeCount = nodesByOrientation.size();
   std::vector<ComplexMatrix> res(3, ComplexMatrix(nodeCount));
   
   for(size_t c = 0; c < 3; c++)
   {
      for(size_t i = 0; i < nodeCount; i++)
      {
         res[c][i].resize(harmonicLength);
      }
   }

   for(size_t h = 0; h < nodeCount; h++)
   {
      auto& expansion = nodesByOrientation[h]->multipoleExpansion();

      for(size_t c = 0; c < 3; c++)
      {
         auto ñomplex = Harmonics::realToComplex(Harmonics::separateCoord(expansion, c));

         for(size_t i = 0; i < harmonicLength; i++)
         {
            res[c][h][i] = ñomplex.getHarmonic(i);
         }
      }
   }

   return res;
}

ComplexMatrix MultipoleSolver::getExpansionsInOneOrientationAsVectors(
   const std::vector<OctreeNode*>& nodesByOrientation,
   size_t padding)
{
   size_t width = math::nextDevisible(
      harmonicLength,
      padding);

   size_t height = math::nextDevisible(
      nodesByOrientation.size(),
      padding);

   ComplexMatrix res(3, std::vector<Complex>(width * height));

   for(int nodeId = 0; nodeId < nodesByOrientation.size(); nodeId++)
   {
      auto& expansion = nodesByOrientation[nodeId]->multipoleExpansion();

      for(size_t c = 0; c < 3; c++)
      {
         auto complex = Harmonics::realToComplex(Harmonics::separateCoord(expansion, c));
         std::copy(complex.data().begin(), 
                   complex.data().end(), 
                   res[c].begin() + nodeId * width);
      }
   }

   return res;
}

void MultipoleSolver::accountChildrenContributions(
   const std::vector<OctreeNode*>& nodesByOrientation,
   const std::vector<ComplexMatrix>& contributions)
{
   for(size_t nodeId = 0; nodeId < nodesByOrientation.size(); nodeId++)
   {
      auto parent = nodesByOrientation[nodeId]->parent();

      auto xComponent = Harmonics::complexToReal(
         ComplexHarmonicSeries(contributions[0][nodeId]));

      auto yComponent = Harmonics::complexToReal(
         ComplexHarmonicSeries(contributions[1][nodeId]));

      auto zComponent = Harmonics::complexToReal(
         ComplexHarmonicSeries(contributions[2][nodeId]));
      
      auto totalContribution = Harmonics::createFormXYZ(
         xComponent, yComponent, zComponent);

      parent->multipoleExpansion().add(totalContribution);
   }
}

void MultipoleSolver::accountChildrenContributions(
   const std::vector<OctreeNode*>& nodesByOrientation,
   const ComplexMatrix& contributions,
   size_t padding)
{
   for(int nodeId = 0; nodeId < nodesByOrientation.size(); nodeId++)
   {
      auto parent = nodesByOrientation[nodeId]->parent();

      auto xComponent = Harmonics::complexToReal(
         math::getRow(
            contributions[0],
            harmonicLength,
            padding,
            nodeId));

      auto yComponent = Harmonics::complexToReal(
         math::getRow(
            contributions[1],
            harmonicLength,
            padding,
            nodeId));

      auto zComponent = Harmonics::complexToReal(
         math::getRow(
            contributions[2],
            harmonicLength,
            padding,
            nodeId));

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
