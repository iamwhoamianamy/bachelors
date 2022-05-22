#include "fast_multipole_solver.hpp"

#include <iostream>

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
      //Box(Vector3(10, 5, 8), Vector3(1, 1, 1)), calculationPointOctreeLeafCapacity);
   _calculationPointOctreeRoot->insert(_points);
}

void FastMultipoleSolver::calcLocalMultipoleExpansionsWithComplexTranslation() const
{
   auto nodesToVisit = _calculationPointOctreeRoot->getAllNodesAsSet();
   auto nodesToAvoid = std::set<CalculationPointOctreeNode*>();

   _quadratureOctreeRoot->translateMultipoleExpansionsToLocal(
      _calculationPointOctreeRoot,
      nodesToVisit,
      nodesToAvoid);

   _calculationPointOctreeRoot->propagateLocalExpansions();
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
      /*HarmonicSeries<Vector3> nearestContribution(harmonicOrder);
      auto nodesInsidePointsBox = _quadratureOctreeRoot->getBiggestNodesInBox(
         node->box());
      
      for (auto node : nodesInsidePointsBox)
      {
         if(2 * node->box().radius() < (point - node->box().center()).length())
         {
            auto irregularHarmonic = Harmonics::calcIrregularSolidHarmonics(
               harmonicOrder, point - node->box().center());

            answer += mult(node->multipoleExpansion(), irregularHarmonic);
         }
      }*/

      /*HarmonicSeries<Vector3> nearestContribution(harmonicOrder);
      auto biggestNearNodes = _quadratureOctreeRoot->getBiggestNearNodes(
         point, node->box().radius());

      for(auto nearNode : biggestNearNodes)
      {
         auto irregularHarmonic = Harmonics::calcIrregularSolidHarmonics(
            harmonicOrder, point - nearNode->box().center());

         answer += mult(nearNode->multipoleExpansion(), irregularHarmonic);
      }*/
      
      /*HarmonicSeries<Vector3> nearestContribution(harmonicOrder);
      auto newLeftToInteractWith = node->leftToInteractWith;

      for(auto nearNode : node->leftToInteractWith)
      {
         real radius1 = nearNode->box().radius();
         real radius2 = node->box().radius();

         real minRadius = std::min(radius1, radius2);
         real maxRadius = std::max(radius1, radius2);

         real minimalDistance = 2 * minRadius + maxRadius;
         real distance1 = (point - nearNode->box().center()).length();
         real distance2 = (node->box().center() - nearNode->box().center()).length();

         if(2 * radius1 < distance1)
         {
            nearNode->removeAllDescendantsFromSet(newLeftToInteractWith);
         }
         else
         {
            newLeftToInteractWith.erase(nearNode);
         }
      }

      for(auto nearNode : newLeftToInteractWith)
      {
         auto irregularHarmonic = Harmonics::calcIrregularSolidHarmonics(
            harmonicOrder, point - nearNode->box().center());

         answer += mult(nearNode->multipoleExpansion(), irregularHarmonic);
      }*/

      result.emplace_back(point, answer / (4.0 * math::PI) * current);
   }
   return result;
}

std::vector<Vector3> FastMultipoleSolver::calcB(real current)
{
   std::vector<Vector3> res(_points.size());



   return res;
}
