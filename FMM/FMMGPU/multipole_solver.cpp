#include <iostream>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <omp.h>

#include "cblas.h"
#include "multipole_solver.hpp"
#include "math.hpp"
#include "integration.hpp"
#include "harmonics.hpp"
#include "math.hpp"
#include "translation_algorithms.hpp"
#include "kernels.cuh"
#include "testing_helpers.hpp"
#include "multipole_translator.hpp"

MultipoleSolver::MultipoleSolver(std::vector<Quadrature>& quadratures,
                                 size_t octreeLeafCapacity) :
   _quadratures(quadratures), octreeLeafCapacity(octreeLeafCapacity)
{
   quadratureOctreeRoot = new QuadratureOctreeNode(Box(Vector3(0, 0, 0), Vector3(3, 3, 3)), octreeLeafCapacity);
   quadratureOctreeRoot->insert(_quadratures);

   _realToComplexMatrix = Harmonics::calcRealToComplexMatrixTransposed1D(
      harmonicOrder);

   _complexToRealMatrix = Harmonics::calcComplexToRealMatrixTransposed1D(
      harmonicOrder);
}

void MultipoleSolver::calcMultipolesAtLeaves()
{
   quadratureOctreeRoot->calcMultipoleExpansionsAtLeaves(harmonicOrder);
   _multipolesAtLeavesAreReady = true;
}

size_t MultipoleSolver::getOctreeNodeCount() const
{
   return quadratureOctreeRoot->getAllNodeCount();
}

void MultipoleSolver::calcLocalMultipoles(M2MAlg algorithm, M2MDevice device)
{
   if(_multipolesAtLeavesAreReady)
   {
      switch(algorithm)
      {
         case M2MAlg::NoTranslation:
         {
            calcLocalMultipolesWithoutTranslation();
            break;
         }
         case M2MAlg::ComplexTranslation:
         {
            calcLocalMultipolesWithComplexTranslation();
            break;
         }
         case M2MAlg::RealTranslation:
         {
            calcLocalMultipolesWithRealTranslation();
            break;
         }
         case M2MAlg::Layers:
         {
            calcLocalMultipolesWithLayersOrMatrices(device, false);
            break;
         }
         case M2MAlg::Matrices:
         {
            calcLocalMultipolesWithLayersOrMatrices(device, true);
            break;
         }
      }
   }
   else
   {
      throw std::exception("Multipoles at leaves are not ready!");
   }
}

void MultipoleSolver::calcLocalMultipolesWithoutTranslation()
{
   quadratureOctreeRoot->calcMultipoleExpansionsWithoutTranslation(harmonicOrder);
   _multipolesAreReady = true;
}

void MultipoleSolver::calcLocalMultipolesWithComplexTranslation()
{
   quadratureOctreeRoot->calcMultipoleExpansionsWithComplexTranslation(harmonicOrder);
   _multipolesAreReady = true;
}

void MultipoleSolver::calcLocalMultipolesWithRealTranslation()
{
   quadratureOctreeRoot->calcMultipoleExpansionsWithRealTranslation(harmonicOrder);
   _multipolesAreReady = true;
}

void MultipoleSolver::calcLocalMultipolesWithLayersOrMatrices(
   M2MDevice device,
   bool useMatrices)
{
   std::vector<std::vector<QuadratureOctreeNode*>> layers;
   enumerateNodes(quadratureOctreeRoot, layers, 0);
   quadratureOctreeRoot->initAllMultipoleExpansions(harmonicOrder);

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
   _multipolesAreReady = true;
}

void MultipoleSolver::enumerateNodes(
   QuadratureOctreeNode* node,
   std::vector<std::vector<QuadratureOctreeNode*>>& layers,
   size_t currentLayerId)
{
   if(!node->quadratures().empty() || node->isSubdivided())
   {
      if(layers.size() <= currentLayerId)
         layers.push_back(std::vector<QuadratureOctreeNode*>());

      layers[currentLayerId].push_back(node);

      for(auto child : node->children())
      {
         enumerateNodes(child, layers, currentLayerId + 1);
      }
   }
}

void MultipoleSolver::calcContributionsToHigherLayers(
   const std::vector<std::vector<QuadratureOctreeNode*>>& layers,
   M2MDevice device,
   bool useMatrices)
{
   useMatrices ? calcContributionsToHigherLevelsWithMatrices(layers, device) :
      calcContributionsToHigherLayers(layers, device);
}

void MultipoleSolver::calcMultipolesAtLeaves(
   const std::vector<std::vector<QuadratureOctreeNode*>>& layers)
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
   const std::vector<std::vector<QuadratureOctreeNode*>>& layers,
   M2MDevice device)
{
   for(int l = layers.size() - 1; l >= 1; l--)
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
   M2MDevice device)
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

   switch(device)
   {
      case M2MDevice::CPU:
      {
         kernels::translateAllCPU(
            result.data(),
            regulars.data(),
            harmonics.data(),
            layer.size(),
            harmonicOrder);

         break;
      }
      case M2MDevice::GPU:
      {
         kernels::translateAllGPU(
            result.data(),
            regulars.data(),
            harmonics.data(),
            layer.size(),
            harmonicOrder);

         break;
      }
      case M2MDevice::Adaptive:
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
   M2MDevice device)
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

            RealMatrix translated(3, std::vector<real>(harmonicLength * nodesCount));

            for(size_t c = 0; c < 3; c++)
            {
               auto kernelStart = std::chrono::steady_clock::now();

               switch(device)
               {
                  case M2MDevice::CPU:
                  {
                     kernels::translateAllCPUMatrixBLAS(
                        translated[c].data(),
                        expansionVectors[c].data(),
                        regularVectors[o].data(),
                        nodesCount,
                        harmonicOrder);

                     break;
                  }
                  case M2MDevice::GPU:
                  {
                    kernels::translateAllGPUMatrixCuBLAS(
                       translated[c].data(),
                        expansionVectors[c].data(),
                        regularVectors[o].data(),
                        nodesCount,
                        harmonicOrder);

                     break;
                  }
                  case M2MDevice::Adaptive:
                  {
                     /*if(expansionVectors[c].size() < adaptiveBorder)
                     {
                        kernels::translateAllCPUMatrixBLAS(
                           t.data(),
                           expansionVectors[c].data(),
                           regularVectors[o].data(),
                           nodesCount,
                           harmonicOrder);
                     }
                     else
                     {
                        kernels::translateAllGPUMatrixCuBLAS(
                           t.data(),
                           expansionVectors[c].data(),
                           regularVectors[o].data(),
                           nodesCount,
                           harmonicOrder);
                     }*/

                     break;
                  }
               }

               auto kernelStop = std::chrono::steady_clock::now();
               kernelTime += std::chrono::duration_cast<std::chrono::microseconds>
                  (kernelStop - kernelStart).count() * 1e-6;
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
            nodesByOrientation[i][0]->box().center - parent->box().center;

         auto regularHarmonics = Harmonics::calcRegularSolidHarmonics(
            harmonicOrder,
            translation);

         res.emplace_back(formMatrixFromRegularHarmonics(
            Harmonics::realToComplex(regularHarmonics)));
      }
   }

   return res;
}

RealMatrix MultipoleSolver::calcRegularMatricesForLayerAsVectors(
   const Matrix<QuadratureOctreeNode*>& nodesByOrientation)
{
   size_t matrixElemCount = harmonicLength * harmonicLength;

   RealMatrix result(8, std::vector<real>(matrixElemCount));

   for(int i = 0; i < 8; i++)
   {
      auto parent = nodesByOrientation[i][0]->parent();

      auto translation = nodesByOrientation[i][0]->box().center -
         parent->box().center;

      auto regularHarmonics = Harmonics::calcRegularSolidHarmonics(
         harmonicOrder, translation);

      auto regularHarmonicsMatrix = formMatrixFromRegularHarmonicsAsVectors(
         Harmonics::realToComplex(regularHarmonics));

      Complex alpha = make_cuComplex(1, 0);
      Complex beta = make_cuComplex(0, 0);

      std::vector<Complex> temp1(matrixElemCount);

      cblas_cgemm(
         CBLAS_ORDER::CblasRowMajor,
         CBLAS_TRANSPOSE::CblasNoTrans,
         CBLAS_TRANSPOSE::CblasNoTrans,
         121, 121, 121,
         (float*)&alpha,
         (float*)_realToComplexMatrix.data(),
         121,
         (float*)regularHarmonicsMatrix.data(),
         121,
         (float*)&beta,
         (float*)temp1.data(),
         121);

      std::vector<Complex> temp2(matrixElemCount);

      cblas_cgemm(
         CBLAS_ORDER::CblasRowMajor,
         CBLAS_TRANSPOSE::CblasNoTrans,
         CBLAS_TRANSPOSE::CblasNoTrans,
         121, 121, 121,
         (float*)&alpha,
         (float*)temp1.data(),
         121,
         (float*)_complexToRealMatrix.data(),
         121,
         (float*)&beta,
         (float*)temp2.data(),
         121);

      cblas_scopy(
         matrixElemCount,
         (float*)(temp2.data()), 2,
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

std::vector<Complex> MultipoleSolver::formMatrixFromRegularHarmonicsAsVectors(
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
                  res[(l * l + l + m) * harmonicLength + (dl * dl + dl + dm)] =
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
         cblas_scopy(
            harmonicLength,
            (float*)expansion.data().data() + c, 3,
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

      std::vector<Vector3> totalContribution(harmonicLength);

      for(size_t c = 0; c < 3; c++)
      {
         cblas_saxpy(
            harmonicLength, 1, 
            contributions[c].data() + harmonicLength * nodeId, 1,
            (float*)(parent->multipoleExpansion().data().data()) + c, 3);
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

Vector3 MultipoleSolver::calcA(real current, const Vector3& point)
{
   if(!_multipolesAreReady)
      throw new std::exception("Multipoles are not ready!");

   return quadratureOctreeRoot->calcA(point) / (4.0 * math::PI) * current;
}

Vector3 MultipoleSolver::calcB(real current, const Vector3& point)
{
   if(!_multipolesAreReady)
      throw new std::exception("Multipoles are not ready!");

   return quadratureOctreeRoot->caclRot(point) / (4.0 * math::PI) * current * math::mu0;
}

MultipoleSolver::~MultipoleSolver()
{
   delete quadratureOctreeRoot;
}