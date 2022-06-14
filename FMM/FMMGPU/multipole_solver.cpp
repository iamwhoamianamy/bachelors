#include <iostream>
#include <chrono>
#include <iomanip>
#include <fstream>

#include "cblas.h"
#include "multipole_solver.hpp"

#include <queue>

#include "math.hpp"
#include "integration.hpp"
#include "harmonics.hpp"
#include "translation_algorithms.hpp"
#include "kernels.cuh"
#include "testing_helpers.hpp"
#include "multipole_translator.hpp"
#include "blass_callers.hpp"

MultipoleSolver::MultipoleSolver(
   std::vector<Quadrature>& quadratures,
   Problem problem,
   size_t quadratureOctreeLeafCapacity) :
   quadratureOctreeLeafCapacity(quadratureOctreeLeafCapacity),
   _problem(problem)
{
   _quadratures.reserve(quadratures.size());

   for(auto& quadrature : quadratures)
   {
      _quadratures.emplace_back(&quadrature);
   }

   initTrees();
   initTransitionMatrices();
}

MultipoleSolver::MultipoleSolver(
   std::vector<BEMQuadrature>& quadratures,
   Problem problem,
   size_t quadratureOctreeLeafCapacity) :
   quadratureOctreeLeafCapacity(quadratureOctreeLeafCapacity),
   _problem(problem)
{
   _quadratures.reserve(quadratures.size());

   for(auto& quadrature : quadratures)
   {
      _quadratures.emplace_back(&quadrature);
   }

   initTrees();
   initTransitionMatrices();
}

void MultipoleSolver::calcMultipoleExpansionsAtLeaves()
{
   std::queue<QuadratureOctreeNode*> qq;

   qq.push(_quadratureOctreeRoot);
   
   while(!qq.empty())
   {
      auto currentNode = qq.front();
      qq.pop();
      
      if(currentNode->isUsefullLeaf())
      {
         switch(_problem)
         {
            case Problem::BioSavartLaplace:
            {
               currentNode->multipoleExpansion() = math::calcIntegralContribution(
                  currentNode->quadratures(), harmonicOrder, currentNode->box().center());
               break;
            }
            case Problem::BEM:
            {
               currentNode->multipoleExpansion() = math::calcBEMIntegralContribution(
                  currentNode->quadratures(), harmonicOrder, currentNode->box().center());
               break;
            }
         }
      }
      else
      {
         for (auto child : currentNode->children())
         {
            qq.push(child);
         }
      }
   }
   
   _multipolesAtLeavesAreReady = true;
}

size_t MultipoleSolver::getQuadratureOctreeNodeCount() const
{
   return _quadratureOctreeRoot->getAllNodeCount();
}

void MultipoleSolver::calclMultipoleExpansions(M2MAlg algorithm, Device device)
{
   if(!_multipolesAtLeavesAreReady)
   {
      calcMultipoleExpansionsAtLeaves();
   }

   switch(algorithm)
   {
      case M2MAlg::NoTranslation:
      {
         calcMultipoleExpansionsWithoutTranslation();
         break;
      }
      case M2MAlg::ComplexTranslation:
      {
         calcMultipoleExpansionsWithComplexTranslation();
         break;
      }
      case M2MAlg::RealTranslation:
      {
         calcMultipoleExpansionsWithRealTranslation();
         break;
      }
      case M2MAlg::Layers:
      {
         calcMultipoleExpanstionsWithLayersOrMatrices(device, false);
         break;
      }
      case M2MAlg::Matrices:
      {
         calcMultipoleExpanstionsWithLayersOrMatrices(device, true);
         break;
      }
   }

   _multipolesAreReady = true;
}

void MultipoleSolver::calcMultipoleExpansionsWithoutTranslation()
{
   std::queue<QuadratureOctreeNode*> qq;

   qq.push(_quadratureOctreeRoot);

   while(!qq.empty())
   {
      auto currentNode = qq.front();
      qq.pop();

      if(!currentNode->isUsefullLeaf())
      {
         currentNode->multipoleExpansion() = math::calcIntegralContribution(
            currentNode->getAllQuadratures(),
            harmonicOrder,
            currentNode->box().center());

         for(auto child : currentNode->children())
         {
            qq.push(child);
         }
      }
   }
}

void MultipoleSolver::calcMultipoleExpansionsWithComplexTranslation()
{
   _quadratureOctreeRoot->calcMultipoleExpansionsWithComplexTranslation(harmonicOrder);
}

void MultipoleSolver::calcMultipoleExpansionsWithRealTranslation()
{
   _quadratureOctreeRoot->calcMultipoleExpansionsWithRealTranslation(harmonicOrder);
}

void MultipoleSolver::calcMultipoleExpanstionsWithLayersOrMatrices(
   Device device,
   bool useMatrices)
{
   std::vector<std::vector<QuadratureOctreeNode*>> layers;
   enumerateNodes(_quadratureOctreeRoot, layers, 0);
   _quadratureOctreeRoot->initAllMultipoleExpansions(harmonicOrder);

   if(log)
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

   calcContributionsToHigherLayers(layers, device, useMatrices);
}

void MultipoleSolver::enumerateNodes(
   QuadratureOctreeNode* node,
   std::vector<std::vector<QuadratureOctreeNode*>>& layers,
   size_t currentLayerId)
{
   if(node->isUsefullLeaf() || node->isSubdivided())
   {
      if(layers.size() <= currentLayerId)
         layers.emplace_back(std::vector<QuadratureOctreeNode*>());

      layers[currentLayerId].push_back(node);

      for(auto child : node->children())
      {
         enumerateNodes(child, layers, currentLayerId + 1);
      }
   }
}

void MultipoleSolver::calcContributionsToHigherLayers(
   const std::vector<std::vector<QuadratureOctreeNode*>>& layers,
   Device device,
   bool useMatrices)
{
   useMatrices ? calcContributionsToHigherLevelsWithMatrices(layers, device) :
      calcContributionsToHigherLayers(layers, device);
}

void MultipoleSolver::calcContributionsToHigherLayers(
   const std::vector<std::vector<QuadratureOctreeNode*>>& layers,
   Device device)
{
   for(size_t l = layers.size() - 1; l >= 1; l--)
   {
      auto start = std::chrono::steady_clock::now();

      if(log)
      {
         std::cout << std::setw(10) << l << std::setw(15) << layers[l].size();
      }

      std::vector<Vector3> contributions =
         calcContributionsToHigherLayer(layers[l], device);

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

      if(log)
      {
         std::cout << std::setw(15) << time << std::endl;
      }
   }
}

std::vector<Vector3> MultipoleSolver::calcContributionsToHigherLayer(
   const std::vector<QuadratureOctreeNode*>& layer,
   Device device)
{
   std::vector<Vector3> harmonics(layer.size() * harmonicLength);
   std::vector<real> regulars(layer.size() * harmonicLength);

   for(size_t i = 0; i < layer.size(); i++)
   {
      Vector3 translation = layer[i]->box().center() - layer[i]->parent()->box().center();
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

   switch(device)
   {
      case Device::CPU:
      {
         kernels::translateAllCPU(
            result.data(),
            regulars.data(),
            harmonics.data(),
            layer.size(),
            harmonicOrder);

         break;
      }
      case Device::GPU:
      {
         kernels::translateAllGPU(
            result.data(),
            regulars.data(),
            harmonics.data(),
            layer.size(),
            harmonicOrder);

         break;
      }
      case Device::Adaptive:
      {
         break;
      }
   }

   auto stop = std::chrono::steady_clock::now();
   double layerTime = test::getTime(start, stop);

   if(log)
   {
      std::cout << std::setw(15) << layerTime;
   }

   return result;
}

void MultipoleSolver::calcContributionsToHigherLevelsWithMatrices(
   const std::vector<std::vector<QuadratureOctreeNode*>>& layers,
   Device device)
{
   for(size_t l = layers.size() - 1; l >= 1; l--)
   {
      double kernelTime = 0;

      auto start = std::chrono::steady_clock::now();

      if(log)
      {
         std::cout << std::setw(10) << l << std::setw(15) << layers[l].size();
      }

      auto& layer = layers[l];
      auto nodesByOrientation = separateNodesByOrientation(layer);
      auto regularVectors = calcRegularMatricesForM2MAsVectors(
         nodesByOrientation);
      
      for(size_t o = 0; o < 8; o++)
      {
         if(!nodesByOrientation[o].empty())
         {
            size_t nodesCount = nodesByOrientation[o].size();

            auto expansionVectors =
               getExpansionsInOneOrientationAsVectors(
                  nodesByOrientation[o]);

            RealMatrix translated(3, std::vector<real>(harmonicLength * nodesCount));

            switch(device)
            {
               case Device::CPU:
               {
                  for(size_t c = 0; c < 3; c++)
                  {
                     auto kernelStart = std::chrono::steady_clock::now();

                     kernels::translateAllCPUMatrixBLAS(
                        translated[c].data(),
                        expansionVectors[c].data(),
                        regularVectors[o].data(),
                        nodesCount,
                        harmonicOrder);

                     auto kernelStop = std::chrono::steady_clock::now();
                     kernelTime += std::chrono::duration_cast<std::chrono::microseconds>
                        (kernelStop - kernelStart).count() * 1e-6;
                  }

                  break;
               }
               case Device::GPU:
               {
                  kernels::translateAllGPUMatrixCuBLAS(
                     translated,
                     expansionVectors,
                     regularVectors[o].data(),
                     nodesCount,
                     harmonicOrder);

                  break;
               }
            }
            
            accountChildrenContributions(
               nodesByOrientation[o],
               translated);
         }
      }

      auto stop = std::chrono::steady_clock::now();
      double layerTime = test::getTime(start, stop);

      if(log)
      {
         std::cout << std::setw(15) << kernelTime;
         std::cout << std::setw(15) << layerTime << std::endl;
      }
   }
}

Matrix<QuadratureOctreeNode*> MultipoleSolver::separateNodesByOrientation(
   const std::vector<QuadratureOctreeNode*>& layer)
{
   Matrix<QuadratureOctreeNode*> res(8);

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
   const Matrix<QuadratureOctreeNode*>& nodesByOrientation)
{
   std::vector<ComplexMatrix> res;
   res.reserve(8);

   for(size_t i = 0; i < 8; i++)
   {
      auto parent = nodesByOrientation[i][0]->parent();

      if(!nodesByOrientation[i].empty())
      {
         auto translation = 
            nodesByOrientation[i][0]->box().center() - parent->box().center();

         auto regularHarmonics = Harmonics::calcRegularSolidHarmonics(
            harmonicOrder,
            translation);

         res.emplace_back(formMatrixFromRegularHarmonics(
            Harmonics::realToComplex(regularHarmonics)));
      }
   }

   return res;
}

RealMatrix MultipoleSolver::calcRegularMatricesForM2MAsVectors(
   const Matrix<QuadratureOctreeNode*>& nodesByOrientation)
{
   size_t matrixElemCount = harmonicLength * harmonicLength;

   RealMatrix result(8, std::vector<real>(matrixElemCount));

   for(int i = 0; i < 8; i++)
   {
      auto parent = nodesByOrientation[i][0]->parent();

      auto translation = nodesByOrientation[i][0]->box().center() -
         parent->box().center();

      auto regularHarmonics = Harmonics::calcRegularSolidHarmonics(
         harmonicOrder, translation);

      auto regularHarmonicsMatrix = formMatrixFromRegularHarmonicsForM2MAsVectors(
         Harmonics::realToComplex(regularHarmonics));

      Complex alpha = makeComplex(1, 0);
      Complex beta = makeComplex(0, 0);

      std::vector<Complex> temp1(matrixElemCount);

      blas::multComplexMatrices(
         CBLAS_ORDER::CblasRowMajor,
         CBLAS_TRANSPOSE::CblasNoTrans,
         CBLAS_TRANSPOSE::CblasNoTrans,
         harmonicLength, harmonicLength, harmonicLength,
         (real*)&alpha,
         (real*)_realToComplexMatrix.data(),
         harmonicLength,
         (real*)regularHarmonicsMatrix.data(),
         harmonicLength,
         (real*)&beta,
         (real*)temp1.data(),
         harmonicLength);

      std::vector<Complex> temp2(matrixElemCount);

      blas::multComplexMatrices(
         CBLAS_ORDER::CblasRowMajor,
         CBLAS_TRANSPOSE::CblasNoTrans,
         CBLAS_TRANSPOSE::CblasNoTrans,
         harmonicLength, harmonicLength, harmonicLength,
         (real*)&alpha,
         (real*)temp1.data(),
         harmonicLength,
         (real*)_complexToRealMatrix.data(),
         harmonicLength,
         (real*)&beta,
         (real*)temp2.data(),
         harmonicLength);

      blas::copyVector(
         matrixElemCount,
         (real*)(temp2.data()), 2,
         result[i].data(), 1);
   }

   return result;
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
                     MultipoleTranslator::multipoleTranslationFactor(m, mu);
               }
            }
         }
      }
   }

   return res;
}

std::vector<Complex> MultipoleSolver::formMatrixFromRegularHarmonicsForM2MAsVectors(
   const ComplexHarmonicSeries& regular)
{
   std::vector<Complex> res(harmonicLength * harmonicLength);

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
                  res[(l * l + l + m) + (dl * dl + dl + dm) * harmonicLength] =
                     regular.getHarmonic(lambda * lambda + lambda + mu) *
                     MultipoleTranslator::multipoleTranslationFactor(m, mu);
               }
            }
         }
      }
   }

   return res;
}

RealMatrix MultipoleSolver::getExpansionsInOneOrientationAsVectors(
   const std::vector<QuadratureOctreeNode*>& nodesByOrientation)
{
   size_t nodeCount = nodesByOrientation.size();

   RealMatrix res(3, std::vector<real>(harmonicLength * nodeCount));

   for(int nodeId = 0; nodeId < nodesByOrientation.size(); nodeId++)
   {
      auto& expansion = nodesByOrientation[nodeId]->multipoleExpansion();

      for(size_t c = 0; c < 3; c++)
      {
         blas::copyVector(
            harmonicLength,
            (real*)expansion.data().data() + c, 3,
            res[c].data() + nodeId * harmonicLength, 1);
      }
   }

   return res;
}

void MultipoleSolver::accountChildrenContributions(
   const std::vector<QuadratureOctreeNode*>& nodesByOrientation,
   const RealMatrix& contributions)
{
   for(int nodeId = 0; nodeId < nodesByOrientation.size(); nodeId++)
   {
      auto parent = nodesByOrientation[nodeId]->parent();

      for(size_t c = 0; c < 3; c++)
      {
         blas::addVectorToVector(
            harmonicLength, 1, 
            contributions[c].data() + harmonicLength * nodeId, 1,
            (real*)(parent->multipoleExpansion().data().data()) + c, 3);
      }
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

void MultipoleSolver::initTrees()
{
   switch(_problem)
   {
      case Problem::BioSavartLaplace:
      {
         _quadratureOctreeRoot = new QuadratureOctreeNode(
            Box(Vector3(0, 0, 0), Vector3(3, 3, 3)), quadratureOctreeLeafCapacity);
         break;
      }
      case Problem::BEM:
      {
         _quadratureOctreeRoot = new QuadratureOctreeNode(
            Box(Vector3(0, 0, 0), Vector3(1.1, 1.1, 1.6)), quadratureOctreeLeafCapacity);
         break;
      }
   }

   _quadratureOctreeRoot->insert(_quadratures);
}

void MultipoleSolver::initTransitionMatrices()
{
   _realToComplexMatrix = Harmonics::calcRealToComplexTransitionMatrix1D(
      harmonicOrder);

   _complexToRealMatrix = Harmonics::calcComplexToRealTransitionMatrix1D(
      harmonicOrder);
}

Vector3 MultipoleSolver::calcA(real current, const Vector3& point)
{
   if(!_multipolesAreReady)
      throw new std::exception("Multipoles are not ready!");

   return _quadratureOctreeRoot->calcA(point) / (4.0 * math::PI) * current;
}

Vector3 MultipoleSolver::calcB(real current, const Vector3& point)
{
   if(!_multipolesAreReady)
      throw new std::exception("Multipoles are not ready!");

   return _quadratureOctreeRoot->calcRot(point) / (4.0 * math::PI) * current * math::MU0;
}

MultipoleSolver::~MultipoleSolver()
{
   delete _quadratureOctreeRoot;
}