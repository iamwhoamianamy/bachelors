#include <iomanip>
#include <iostream>
#include <chrono>

#include "fast_multipole_solver.hpp"
#include "multipole_translator.hpp"
#include <cblas.h>

#include "testing_helpers.hpp"
#include "translation_algorithms.hpp"

FastMultipoleSolver::FastMultipoleSolver(
   std::vector<Quadrature>& quadratures,
   std::vector<Vector3>& points,
   size_t quadratureOctreeLeafCapacity,
   size_t calculationPointOctreeLeafCapacity) :
   MultipoleSolver(quadratures, quadratureOctreeLeafCapacity),
   _points(points),
   calculationPointOctreeLeafCapacity(calculationPointOctreeLeafCapacity)
{
   initTrees();
}

FastMultipoleSolver::~FastMultipoleSolver()
{
   delete _calculationPointOctreeRoot;
}

void FastMultipoleSolver::initTrees()
{
   MultipoleSolver::initTrees();
   _calculationPointOctreeRoot = new CalculationPointOctreeNode(
      math::getBoundingBox(_points), calculationPointOctreeLeafCapacity);
   _calculationPointOctreeRoot->insert(_points);
}

void FastMultipoleSolver::calclLocalMultipoleExpansions(
   M2LAlg algorithm,
   Device device)
{
   if(_multipolesAreReady)
   {
      if(!_localMultipolesAreInitialized)
      {
         _calculationPointOctreeRoot->initAllLocalExpansions(harmonicOrder);
         _localMultipolesAreInitialized = true;
      }

      formInteractionMaps(
         _quadratureOctreeRoot,
         _calculationPointOctreeRoot);

      accountFarInteractions();

      switch(algorithm)
      {
         case M2LAlg::ComplexTranslation:
         {
            propagateLocalExpansionsWithComplexTranslation();
            break;
         }
         case M2LAlg::Matrices:
         {
            propagateLocalExpansionsWithMatrices(device);
            break;
         }
         default:
         {
            throw std::exception("Not implemented!");
         }
      }
   }
   else
   {
      throw std::exception("Multipoles are not ready!");
   }
}


void FastMultipoleSolver::formInteractionMaps(
   QuadratureOctreeNode* quadratureTreeNode,
   CalculationPointOctreeNode* calcPointTreeNode)
{
   if(!((quadratureTreeNode->isSubdivided() || quadratureTreeNode->isUsefullLeaf()) &&
        (calcPointTreeNode->isSubdivided() || calcPointTreeNode->isUsefullLeaf())))
   {
      return;
   }

   if(checkIfFarEnough(quadratureTreeNode, calcPointTreeNode))
   {
      _farInteractionMap[calcPointTreeNode].insert(quadratureTreeNode);
   }
   else
   {
      if(calcPointTreeNode->isUsefullLeaf())
      {
         bool сhildrenAreCloseEnough = true;

         for(auto child : quadratureTreeNode->children())
         {
            if(checkIfFarEnough(child, calcPointTreeNode))
            {
               сhildrenAreCloseEnough = false;
               break;
            }
         }

         if(сhildrenAreCloseEnough)
         {
            _closeInteractionMap[calcPointTreeNode].insert(quadratureTreeNode);
         }
         else
         {
            for(auto child : quadratureTreeNode->children())
            {
               formInteractionMaps(child, calcPointTreeNode);
            }
         }
      }
      else
      {
         if((calcPointTreeNode->box().radius() < quadratureTreeNode->box().radius() ||
             calcPointTreeNode->isUsefullLeaf()) && !quadratureTreeNode->isUsefullLeaf())
         {
            for(auto child : quadratureTreeNode->children())
            {
               formInteractionMaps(child, calcPointTreeNode);
            }
         }
         else
         {
            for(auto child : calcPointTreeNode->children())
            {
               formInteractionMaps(quadratureTreeNode, child);
            }
         }
      }
   }
}

bool FastMultipoleSolver::checkIfFarEnough(
   const QuadratureOctreeNode* quadratureTreeNode,
   const CalculationPointOctreeNode* calcPointTreeNode) const
{
   real radius1 = quadratureTreeNode->box().radius();
   real radius2 = calcPointTreeNode->box().radius();

   real minRadius = std::min(radius1, radius2);
   real maxRadius = std::max(radius1, radius2);

   real minimalDistance = 2 * maxRadius + minRadius;

   real distanceSquared = (
      quadratureTreeNode->box().center() -
      calcPointTreeNode->box().center()).lengthSquared();

   return minimalDistance * minimalDistance < distanceSquared;
}

Vector3 FastMultipoleSolver::calcA(real current, const Vector3& point)
{
   return MultipoleSolver::calcA(current, point);
}

Vector3 FastMultipoleSolver::calcB(real current, const Vector3& point)
{
   return MultipoleSolver::calcB(current, point);
}

std::vector<std::pair<Vector3, Vector3>> FastMultipoleSolver::calcA(real current)
{
   auto ffmResults = _calculationPointOctreeRoot->calcA(_points.size());
   std::vector<std::pair<Vector3, Vector3>> result;
   result.reserve(ffmResults.size());

   for(auto& [point, answer, node] : ffmResults)
   {
      for(auto interactionNode : _closeInteractionMap[node])
      {
         answer += interactionNode->calcA(point);
      }

      result.emplace_back(point, answer / (4.0 * math::PI) * current);
   }

   return result;
}

std::vector<std::pair<Vector3, Vector3>> FastMultipoleSolver::calcB(real current)
{
   auto ffmResults = _calculationPointOctreeRoot->calcRot(_points.size());
   std::vector<std::pair<Vector3, Vector3>> result;
   result.reserve(ffmResults.size());

   for(auto& [point, answer, node] : ffmResults)
   {
      for(auto interactionNode : _closeInteractionMap[node])
      {
         answer += interactionNode->calcRot(point);
      }

      result.emplace_back(point, answer / (4.0 * math::PI) * current * math::MU0);
   }

   return result;
}


void FastMultipoleSolver::propagateLocalExpansionsWithComplexTranslation()
{
   _calculationPointOctreeRoot->propagateLocalExpansions();
}

void FastMultipoleSolver::accountFarInteractions()
{
   for(auto& [calcTreeNode, interactionSet] : _farInteractionMap)
   {
      for(auto interactionNode : interactionSet)
      {
         auto translation = interactionNode->box().center() -
            calcTreeNode->box().center();

         auto contributionToLocalExpansion =
            MultipoleTranslator::multipoleToLocalWithComplex(
               interactionNode->multipoleExpansion(),
               translation);

         calcTreeNode->localExpansion().add(
            contributionToLocalExpansion);
      }
   }
}

void FastMultipoleSolver::propagateLocalExpansionsWithMatrices(Device device)
{
   std::vector<std::vector<CalculationPointOctreeNode*>> layers;
   enumerateNodes(_calculationPointOctreeRoot, layers, 0);

   if(log)
   {
      std::cout << std::setw(10) << "layer" << std::setw(15) << "mlpl count";
      std::cout << std::setw(15) << "kernel time" << std::setw(15) << "total time" << std::endl;
      std::cout << "-------------------------------------------------------" << std::endl;
      std::cout << std::fixed;
   }

   calcContributionsToLowerLayers(layers, device);
   //_multipolesAreReady = true;
}

void FastMultipoleSolver::calcContributionsToLowerLayers(
   const std::vector<std::vector<CalculationPointOctreeNode*>>& layers,
   Device device)
{
   for(size_t l = 1; l < layers.size(); l++)
   {
      double kernelTime = 0;

      auto start = std::chrono::steady_clock::now();

      if(log)
      {
         std::cout << std::setw(10) << l << std::setw(15) << layers[l].size();
      }

      auto& layer = layers[l];
      auto nodesByOrientation = separateNodesByOrientation(layer);
      auto regularVectors = calcRegularMatricesForL2LAsVectors(
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
                  case Device::CPU:
                  {
                     kernels::translateAllCPUMatrixBLAS(
                        translated[c].data(),
                        expansionVectors[c].data(),
                        regularVectors[o].data(),
                        nodesCount,
                        harmonicOrder);

                     break;
                  }
                  case Device::GPU:
                  {
                     kernels::translateAllGPUMatrixCuBLAS(
                        translated[c].data(),
                        expansionVectors[c].data(),
                        regularVectors[o].data(),
                        nodesCount,
                        harmonicOrder);

                     break;
                  }
                  default:
                  {
                     throw std::exception("Not implemented!");
                  }
               }

               auto kernelStop = std::chrono::steady_clock::now();
               kernelTime += std::chrono::duration_cast<std::chrono::microseconds>
                  (kernelStop - kernelStart).count() * 1e-6;
            }

            accountContributionsToChildren(
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

void FastMultipoleSolver::accountContributionsToChildren(
   const std::vector<CalculationPointOctreeNode*>& nodesByOrientation,
   const RealMatrix& contributions)
{
   for(int nodeId = 0; nodeId < nodesByOrientation.size(); nodeId++)
   {
      for(size_t c = 0; c < 3; c++)
      {
         cblas_saxpy(
            harmonicLength, 1,
            contributions[c].data() + harmonicLength * nodeId, 1,
            (float*)(nodesByOrientation[nodeId]->localExpansion().data().data()) + c, 3);
      }
   }
}

RealMatrix FastMultipoleSolver::getExpansionsInOneOrientationAsVectors(
   const std::vector<CalculationPointOctreeNode*>& nodesByOrientation) const
{
   size_t nodeCount = nodesByOrientation.size();

   RealMatrix res(3, std::vector<real>(harmonicLength * nodeCount));

   for(int nodeId = 0; nodeId < nodesByOrientation.size(); nodeId++)
   {
      auto& expansion = nodesByOrientation[nodeId]->parent()->localExpansion();

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

void FastMultipoleSolver::enumerateNodes(
   QuadratureOctreeNode* node,
   std::vector<std::vector<QuadratureOctreeNode*>>& layers,
   size_t currentLayerId)
{
   MultipoleSolver::enumerateNodes(node, layers, currentLayerId);
}

void FastMultipoleSolver::enumerateNodes(
   CalculationPointOctreeNode* node,
   std::vector<std::vector<CalculationPointOctreeNode*>>& layers,
   size_t currentLayerId)
{
   if(node->isSubdivided() || node->isUsefullLeaf())
   {
      if(layers.size() <= currentLayerId)
         layers.push_back(std::vector<CalculationPointOctreeNode*>());

      layers[currentLayerId].push_back(node);

      for(auto child : node->children())
      {
         enumerateNodes(child, layers, currentLayerId + 1);
      }
   }
}

Matrix<QuadratureOctreeNode*> FastMultipoleSolver::separateNodesByOrientation(
   const std::vector<QuadratureOctreeNode*>& layer)
{
   return MultipoleSolver::separateNodesByOrientation(layer);
}

Matrix<CalculationPointOctreeNode*> FastMultipoleSolver::separateNodesByOrientation(
   const std::vector<CalculationPointOctreeNode*>& layer)
{
   Matrix<CalculationPointOctreeNode*> res(8);

   /*for(auto node : layer)
   {
      if(node->isSubdivided())
      {
         for(size_t i = 0; i < 8; i++)
         {
            auto child = node->children()[i];

            if(child->isSubdivided() ||
               child->isUsefullLeaf())
            {
               res[i].push_back(node);
            }
         }
      }
   }*/

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

RealMatrix FastMultipoleSolver::calcRegularMatricesForL2LAsVectors(
   const Matrix<CalculationPointOctreeNode*>& nodesByOrientation)
{
   size_t matrixElemCount = harmonicLength * harmonicLength;

   RealMatrix result(8, std::vector<real>(matrixElemCount));

   for(int i = 0; i < 8; i++)
   {
      auto parent = nodesByOrientation[i][0]->parent();

      auto translation = parent->box().center() - 
         nodesByOrientation[i][0]->box().center();

      auto regularHarmonics = Harmonics::calcRegularSolidHarmonics(
         harmonicOrder, translation);

      auto regularHarmonicsMatrix = formMatrixFromRegularHarmonicsForL2LAsVectors(
         Harmonics::realToComplex(regularHarmonics));

      Complex alpha = make_cuComplex(1, 0);
      Complex beta = make_cuComplex(0, 0);

      std::vector<Complex> temp1(matrixElemCount);

      cblas_cgemm(
         CBLAS_ORDER::CblasRowMajor,
         CBLAS_TRANSPOSE::CblasNoTrans,
         CBLAS_TRANSPOSE::CblasNoTrans,
         harmonicLength, harmonicLength, harmonicLength,
         (float*)&alpha,
         (float*)_realToComplexMatrix.data(),
         harmonicLength,
         (float*)regularHarmonicsMatrix.data(),
         harmonicLength,
         (float*)&beta,
         (float*)temp1.data(),
         harmonicLength);

      std::vector<Complex> temp2(matrixElemCount);

      cblas_cgemm(
         CBLAS_ORDER::CblasRowMajor,
         CBLAS_TRANSPOSE::CblasNoTrans,
         CBLAS_TRANSPOSE::CblasNoTrans,
         harmonicLength, harmonicLength, harmonicLength,
         (float*)&alpha,
         (float*)temp1.data(),
         harmonicLength,
         (float*)_complexToRealMatrix.data(),
         harmonicLength,
         (float*)&beta,
         (float*)temp2.data(),
         harmonicLength);

      cblas_scopy(
         matrixElemCount,
         (float*)(temp2.data()), 2,
         result[i].data(), 1);
   }

   return result;
}

std::vector<Complex> FastMultipoleSolver::formMatrixFromRegularHarmonicsForL2LAsVectors(
   const ComplexHarmonicSeries& regular) const
{
   std::vector<Complex> res(harmonicLength * harmonicLength);

   for(int l = 0; l <= regular.order(); l++)
   {
      for(int m = -l; m <= l; m++)
      {
         for(int lambda = 0; lambda <= l; lambda++)
         {
            int dl = lambda - l;

            for(int mu = -lambda; mu <= lambda; mu++)
            {
               int dm = m - mu;

               if(-dl <= dm && dm <= dl && dl <= harmonicOrder)
               {
                  res[(l * l + l + m) + (lambda * lambda + lambda + mu) * harmonicLength] =
                     regular.getHarmonic(dl * dl + dl + dm) *
                     MultipoleTranslator::localTranslationFactor(m, mu, lambda, l);
               }
            }
         }
      }
   }

   return res;
}
