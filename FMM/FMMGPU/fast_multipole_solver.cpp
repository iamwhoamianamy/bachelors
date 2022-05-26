#include "fast_multipole_solver.hpp"

#include <iostream>

#include "multipole_translator.hpp"

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

   real maxX = math::max(_points, 0);
   real maxY = math::max(_points, 1);
   real maxZ = math::max(_points, 2);

   real minX = math::min(_points, 0);
   real minY = math::min(_points, 1);
   real minZ = math::min(_points, 2);

   Vector3 halfDim(
      (maxX - minX) / 2 + 1e-5,
      (maxY - minY) / 2 + 1e-5,
      (maxZ - minZ) / 2 + 1e-5);

   Vector3 center(minX + halfDim.x, minY + halfDim.y, minZ + halfDim.z);

   _calculationPointOctreeRoot = new CalculationPointOctreeNode(
      Box(center, halfDim * 1.1), calculationPointOctreeLeafCapacity);
   _calculationPointOctreeRoot->insert(_points);
}

void FastMultipoleSolver::calcLocalMultipoleExpansionsWithComplexTranslation()
{
   formInteractionMaps(
      _quadratureOctreeRoot,
      _calculationPointOctreeRoot);

   accountFarInteractions();

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

         for (auto child : quadratureTreeNode->children())
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

void FastMultipoleSolver::calclLocalMultipoleExpansions(
   M2LAlg algorithm, 
   M2MDevice device)
{
   if(_multipolesAreReady)
   {
      if(!_localMultipolesAreInitialized)
      {
         _calculationPointOctreeRoot->initAllLocalExpansions(harmonicOrder);
         _localMultipolesAreInitialized = true;
      }

      switch(algorithm)
      {
         case M2LAlg::ComplexTranslation:
         {
            calcLocalMultipoleExpansionsWithComplexTranslation();
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
   auto M2MAndM2LResults = _calculationPointOctreeRoot->calcA(_points.size());
   std::vector<std::pair<Vector3, Vector3>> result;
   result.reserve(M2MAndM2LResults.size());

   for (auto &[point, answer, node] : M2MAndM2LResults)
   {
      for (auto interactionNode : _closeInteractionMap[node])
      {
         answer += interactionNode->calcA(point);
      }

      result.emplace_back(point, answer / (4.0 * math::PI) * current);
   }

   return result;
}

std::vector<std::pair<Vector3, Vector3>> FastMultipoleSolver::calcB(real current)
{
   auto M2MAndM2LResults = _calculationPointOctreeRoot->calcRot(_points.size());
   std::vector<std::pair<Vector3, Vector3>> result;
   result.reserve(M2MAndM2LResults.size());

   for(auto& [point, answer, node] : M2MAndM2LResults)
   {
      for(auto interactionNode : _closeInteractionMap[node])
      {
         answer += interactionNode->caclRot(point);
      }

      result.emplace_back(point, answer / (4.0 * math::PI) * current * math::MU0);
   }

   return result;
}
