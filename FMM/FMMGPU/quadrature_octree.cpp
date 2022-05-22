#include "quadrature_octree.hpp"

#include <iostream>

#include "math.hpp"
#include "integration.hpp"
#include "harmonics.hpp"
#include "multipole_translator.hpp"

#include <thrust/complex.h>
#include <queue>

QuadratureOctreeNode::QuadratureOctreeNode() : 
   _parent(nullptr)
{
}

QuadratureOctreeNode::QuadratureOctreeNode(
   const Box& box,
   const size_t capacity,
   QuadratureOctreeNode* parent) :
   _parent(parent), _box(box), _capacity(capacity)
{

}

void QuadratureOctreeNode::insert(std::vector<Quadrature>& points)
{
   for(auto& point : points)
   {
      insert(point);
   }
}

void QuadratureOctreeNode::insert(Quadrature& point)
{
   if(!_box.contains(point.coordinates))
   {
      return;
   }

   if(_quadratures.size() < _capacity)
   {
      if(!isSubdivided())
      {
         _quadratures.push_back(&point);
      }
      else
      {
         for(auto& child : _children)
         {
            child->insert(point);
         }
      }
   }
   else
   {
      subdivide();

      _quadratures.push_back(&point);

      for(auto p : _quadratures)
      {
         for(auto& child : _children)
         {
            child->insert(*p);
         }
      }

      _quadratures.resize(0);
   }
}

void QuadratureOctreeNode::subdivide()
{
   float x = _box.center().x;
   float y = _box.center().y;
   float z = _box.center().z;

   float w = _box.halfDimensions().x;
   float h = _box.halfDimensions().y;
   float d = _box.halfDimensions().z;

   Vector3 childrenHalfDimensions = { w / 2, h / 2, d / 2 };

   _children.reserve(8);

   _children.emplace_back(new QuadratureOctreeNode(Box({ x - w / 2, y - h / 2, z - d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new QuadratureOctreeNode(Box({ x + w / 2, y - h / 2, z - d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new QuadratureOctreeNode(Box({ x + w / 2, y + h / 2, z - d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new QuadratureOctreeNode(Box({ x - w / 2, y + h / 2, z - d / 2 }, childrenHalfDimensions), _capacity, this));

   _children.emplace_back(new QuadratureOctreeNode(Box({ x - w / 2, y - h / 2, z + d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new QuadratureOctreeNode(Box({ x + w / 2, y - h / 2, z + d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new QuadratureOctreeNode(Box({ x + w / 2, y + h / 2, z + d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new QuadratureOctreeNode(Box({ x - w / 2, y + h / 2, z + d / 2 }, childrenHalfDimensions), _capacity, this));
}

std::vector<Quadrature*> QuadratureOctreeNode::getAllQuadratures() const
{
   std::vector<Quadrature*> result;
   result.insert(result.end(), _quadratures.begin(), _quadratures.end());

   if(isSubdivided())
   {
      for(auto child : _children)
      {
         auto childQudratures = child->getAllQuadratures();
         result.insert(result.end(), childQudratures.begin(), childQudratures.end());
      }
   }

   return result;
}

void QuadratureOctreeNode::calcMultipoleExpansionsWithoutTranslation(int n)
{
   if(!isLeafAndUseful())
   {
      _multipoleExpansion = math::calcIntegralContribution(getAllQuadratures(), n, _box.center());

      for(auto& child : _children)
      {
         child->calcMultipoleExpansionsWithoutTranslation(n);
      }
   }
}

void QuadratureOctreeNode::calcMultipoleExpansionsWithComplexTranslation(int n)
{
   if(!isLeafAndUseful())
   {
      _multipoleExpansion = HarmonicSeries<Vector3>(n);

      for(auto child : _children)
      {
         child->calcMultipoleExpansionsWithComplexTranslation(n);
      }

      for(auto child : _children)
      {
         if(!child->_quadratures.empty() && !child->isSubdivided() ||
            child->_quadratures.empty() && child->isSubdivided())
            _multipoleExpansion.add(MultipoleTranslator::translateMultipoleWithComplex(
               child->multipoleExpansion(), child->box().center() - _box.center()));
      }
   }
}

void QuadratureOctreeNode::calcMultipoleExpansionsWithRealTranslation(int n)
{
   if(!isLeafAndUseful())
   {
      _multipoleExpansion = HarmonicSeries<Vector3>(n);

      for(auto child : _children)
      {
         child->calcMultipoleExpansionsWithRealTranslation(n);
      }

      for(auto child : _children)
      {
         if(!child->_quadratures.empty() && !child->isSubdivided() ||
            child->_quadratures.empty() && child->isSubdivided())
            _multipoleExpansion.add(MultipoleTranslator::translateMultipoleWithReal(
               child->multipoleExpansion(), child->box().center() - _box.center()));
      }
   }
}

void QuadratureOctreeNode::initAllMultipoleExpansions(size_t n)
{
   if(isSubdivided() || _quadratures.empty())
   {
      _multipoleExpansion = HarmonicSeries<Vector3>(n);

      for(auto child : _children)
      {
         child->initAllMultipoleExpansions(n);
      }
   }
}

void QuadratureOctreeNode::calcMultipoleExpansionsAtLeaves(size_t n)
{
   if(isLeafAndUseful())
   {
      _multipoleExpansion = math::calcIntegralContribution(
         _quadratures, n, _box.center());
   }
   else
   {
      for(auto child : _children)
      {
         child->calcMultipoleExpansionsAtLeaves(n);
      }
   }
}

Vector3 QuadratureOctreeNode::calcA(const Vector3& point) const
{
   int n = _multipoleExpansion.order();

   Vector3 res;

   if(2 * _box.radius() < (point - _box.center()).length())
   {
      auto irregularHarmonic = Harmonics::calcIrregularSolidHarmonics(n, point - _box.center());

      for(int l = 0; l < n; l++)
      {
         for(int m = -l; m <= l; m++)
         {
            res += _multipoleExpansion.getHarmonic(l, m) * irregularHarmonic.getHarmonic(l, m);
         }
      }
   }
   else
   {
      for(auto child : _children)
      {
         res += child->calcA(point);
      }
   }

   return res;
}

Vector3 QuadratureOctreeNode::caclRot(const Vector3& point) const
{
   int n = _multipoleExpansion.order();
   real eps = 1e-3;

   Vector3 res;

   if(2 * _box.radius() < (point - _box.center()).length())
   {
      auto hx1 = Harmonics::calcIrregularSolidHarmonics(n, point - _box.center() +
                                                        Vector3::xAxis() * eps);
      auto hx2 = Harmonics::calcIrregularSolidHarmonics(n, point - _box.center() -
                                                        Vector3::xAxis() * eps);

      auto hy1 = Harmonics::calcIrregularSolidHarmonics(n, point - _box.center() +
                                                        Vector3::yAxis() * eps);
      auto hy2 = Harmonics::calcIrregularSolidHarmonics(n, point - _box.center() -
                                                        Vector3::yAxis() * eps);

      auto hz1 = Harmonics::calcIrregularSolidHarmonics(n, point - _box.center() +
                                                        Vector3::zAxis() * eps);
      auto hz2 = Harmonics::calcIrregularSolidHarmonics(n, point - _box.center() -
                                                        Vector3::zAxis() * eps);

      hx1.subtract(hx2);
      hy1.subtract(hy2);
      hz1.subtract(hz2);

      Vector3 tempRes;

      for(int l = 0; l < n; l++)
      {
         for(int m = -l; m <= l; m++)
         {
            tempRes.x += _multipoleExpansion.getHarmonic(l, m).z * hy1.getHarmonic(l, m) -
                         _multipoleExpansion.getHarmonic(l, m).y * hz1.getHarmonic(l, m);

            tempRes.y += _multipoleExpansion.getHarmonic(l, m).x * hz1.getHarmonic(l, m) -
                         _multipoleExpansion.getHarmonic(l, m).z * hx1.getHarmonic(l, m);

            tempRes.z += _multipoleExpansion.getHarmonic(l, m).y * hx1.getHarmonic(l, m) -
                         _multipoleExpansion.getHarmonic(l, m).x * hy1.getHarmonic(l, m);
         }
      }

      res += tempRes / (2 * eps);
   }
   else
   {
      for(auto child : _children)
      {
         res += child->caclRot(point);
      }
   }

   return res;
}

size_t QuadratureOctreeNode::getAllNodeCount() const
{
   size_t count = 1;

   if(isSubdivided())
   {
      count += 8;

      for(auto child : _children)
      {
         count += child->getAllNodeCount();
      }
   }
   
   return count;
}

std::set<QuadratureOctreeNode*> QuadratureOctreeNode::getAllUsefulNodes()
{
   std::set<QuadratureOctreeNode*>res;
   getAllUsefulNodes(res);
   return res;
}

void QuadratureOctreeNode::getAllUsefulNodes(
   std::set<QuadratureOctreeNode*>& res)
{
   if(isSubdivided() || !_quadratures.empty())
   {
      res.insert(this);

      for (auto child : _children)
      {
         child->getAllUsefulNodes(res);
      }
   }
}

void QuadratureOctreeNode::translateMultipoleExpansionsToLocal(
   CalculationPointOctreeNode* calculationPointOctreeRoot,
   std::set<CalculationPointOctreeNode*>& nodesToVisit,
   std::set<CalculationPointOctreeNode*>& nodesToAvoid)
{
   auto interactionList = getInteractionList(calculationPointOctreeRoot, nodesToAvoid);

   if(!interactionList.empty())
   {
      for(auto interactionNode : interactionList)
      {
         auto foundToVisit = nodesToVisit.find(interactionNode);

         if(foundToVisit != nodesToVisit.end())
         {
            nodesToVisit.erase(foundToVisit);
            interactionNode->removeAllDescendantsFromSet(nodesToVisit);

            std::set<CalculationPointOctreeNode*> nodesToErase;
            interactionNode->addMeAndAllParentsToSet(nodesToErase);

            for (auto toErase : nodesToErase)
            {
               auto foundToErase = nodesToVisit.find(toErase);
               if(foundToErase != nodesToVisit.end())
                  nodesToVisit.erase(foundToErase);
            }

            nodesToAvoid.insert(nodesToErase.begin(), nodesToErase.end());

            auto translation = _box.center() - interactionNode->box().center();
            auto contributionToLocalExpansion =
               MultipoleTranslator::multipoleToLocalWithComplex(
                  _multipoleExpansion,
                  translation);
            interactionNode->localExpansion().add(
               contributionToLocalExpansion);
         }
      }
   }

   if(!nodesToVisit.empty())
   {
      for(auto child : _children)
      {
         if(child->isSubdivided() || !child->quadratures().empty())
         {
            std::set<CalculationPointOctreeNode*> newNodesToVisit(nodesToVisit);
            std::set<CalculationPointOctreeNode*> newNodesToAvoid(nodesToAvoid);

            child->translateMultipoleExpansionsToLocal(
               calculationPointOctreeRoot,
               newNodesToVisit,
               newNodesToAvoid);
         }
      }
   }

   for (auto node : nodesToVisit)
   {
      if(!node->points().empty())
         node->leftToInteractWith.insert(this);
   }
}

std::vector<CalculationPointOctreeNode*> QuadratureOctreeNode::getInteractionList(
   CalculationPointOctreeNode* calculationPointOctreeNode,
   const std::set<CalculationPointOctreeNode*>& nodesToAvoid) const
{
   std::vector<CalculationPointOctreeNode*> res;
   std::queue<CalculationPointOctreeNode*> queue;
   queue.push(calculationPointOctreeNode);

   while(!queue.empty())
   {
      auto calcPointTreeNode = queue.front();
      queue.pop();

      if(!nodesToAvoid.contains(calcPointTreeNode))
      {
         real radius1 = box().radius();
         real radius2 = calcPointTreeNode->box().radius();

         real minRadius = std::min(radius1, radius2);
         real maxRadius = std::max(radius1, radius2);

         real minimalDistance = 2 * minRadius + maxRadius;
         real distance = sqrt(Vector3::distanceSquared(
            calcPointTreeNode->box().center(),
            _box.center()));
      
         if(minimalDistance < distance)
         {
            res.push_back(calcPointTreeNode);
         }
         else
         {
            for(auto child : calcPointTreeNode->children())
            {
               if(child->isSubdivided() || !child->points().empty())
               {
                  queue.push(child);
               }
            }
         }
      }
      else
      {
         for(auto child : calcPointTreeNode->children())
         {
            if(child->isSubdivided() || !child->points().empty())
            {
               queue.push(child);
            }
         }
      }
   }

   return res;
}

std::vector<QuadratureOctreeNode*> QuadratureOctreeNode::getBiggestNodesInBox(
   const Box& box)
{
   std::vector<QuadratureOctreeNode*> res;
   getBiggestNodesInBox(box, res);
   return res;
}

std::vector<QuadratureOctreeNode*> QuadratureOctreeNode::getBiggestNearNodes(
   const Vector3& point,
   real radius)
{
   std::vector<QuadratureOctreeNode*> res;
   getBiggestNearNodes(res, point, radius);

   return res;
}

void QuadratureOctreeNode::removeAllDescendantsFromSet(
   std::set<QuadratureOctreeNode*>& set) const
{
   for(auto child : _children)
   {
      set.erase(child);
      child->removeAllDescendantsFromSet(set);
   }
}

void QuadratureOctreeNode::getBiggestNearNodes(
   std::vector<QuadratureOctreeNode*>& res,
   const Vector3& point,
   real radius)
{
   real distance = (point - box().center()).length();

   if(distance < 2 * radius && 2 * box().radius() < distance)
   {
      res.push_back(this);
   }
   else
   {
      for(auto child : _children)
      {
         if(child->isSubdivided() || !child->quadratures().empty())
         {
            child->getBiggestNearNodes(res, point, radius);
         }
      }
   }
}

std::vector<QuadratureOctreeNode*> QuadratureOctreeNode::getNearNeighbours()
{
   std::vector<QuadratureOctreeNode*> res;

   if(parent())
   {
      auto grandparent = parent()->parent();

      if(grandparent)
      {
         auto &parentNearNeighbours = grandparent->children();

         for(auto parentNearNeighbour : parentNearNeighbours)
         {
            for(auto parentNearNeighbourChild : parentNearNeighbour->children())
            {
               res.push_back(parentNearNeighbourChild);
            }
         }
      }
   }

   return res;
}

void QuadratureOctreeNode::getBiggestNodesInBox(
   const Box& box, 
   std::vector<QuadratureOctreeNode*>& res)
{
   if((isSubdivided() || !_quadratures.empty()) &&
      (box.contains(_box) || box.intersects(_box)))
   {
      res.push_back(this);
   }
   else
   {
      for (auto child : _children)
      {
         child->getBiggestNodesInBox(box, res);

      }
   }
}

const Box& QuadratureOctreeNode::box() const
{
   return this->_box;
}

bool QuadratureOctreeNode::isSubdivided() const
{
   return !_children.empty();
}

QuadratureOctreeNode* QuadratureOctreeNode::parent()
{
   return _parent;
}

QuadratureOctreeNode* QuadratureOctreeNode::parent() const
{
   return _parent;
}

std::vector<QuadratureOctreeNode*>& QuadratureOctreeNode::children()
{
   return _children;
}

const std::vector<QuadratureOctreeNode*>& QuadratureOctreeNode::children() const
{
   return _children;
}

const std::vector<Quadrature*>& QuadratureOctreeNode::quadratures() const
{
   return _quadratures;
}

HarmonicSeries<Vector3>& QuadratureOctreeNode::multipoleExpansion()
{
   return _multipoleExpansion;
}

const HarmonicSeries<Vector3>& QuadratureOctreeNode::multipoleExpansion() const
{
   return _multipoleExpansion;
}

QuadratureOctreeNode::~QuadratureOctreeNode()
{
   if(isSubdivided())
   {
      for(auto child : _children)
      {
         delete child;
      }
   }
}

bool QuadratureOctreeNode::isLeafAndUseful() const
{
   return !isSubdivided() && !_quadratures.empty();
}
