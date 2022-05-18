#include "quadrature_octree.hpp"
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

      _quadratures.clear();
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

void QuadratureOctreeNode::quarry(const Box& range, std::vector<Quadrature*>& found)
{
   if(!this->_box.intersects(range))
   {
      return;
   }
   else
   {
      for(auto point : _quadratures)
      {
         if(range.contains(point->coordinates))
         {
            found.push_back(point);
         }
      }

      if(isSubdivided())
      {
         for(auto& child : _children)
         {
            child->quarry(range, found);
         }
      }
   }
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

void QuadratureOctreeNode::translateMultipoleExpansionsToLocal(
   CalculationPointOctreeNode* calculationPointOctreeRoot,
   std::set<CalculationPointOctreeNode*>& nodesToVisit)
{
   auto interactionList = getInteractionList(calculationPointOctreeRoot);

   if(!interactionList.empty())
   {
      for(auto interactionNode : interactionList)
      {
         auto found = nodesToVisit.find(interactionNode);

         if(found != nodesToVisit.end())
         {
            nodesToVisit.erase(found);

            auto translation = _box.center() - interactionNode->box().center();
            interactionNode->localExpansion().add(
               MultipoleTranslator::multipoleToLocalWithComplex(
                  _multipoleExpansion,
                  translation));
         }
      }

      if(!nodesToVisit.empty())
      {
         for(auto child : _children)
         {
            child->translateMultipoleExpansionsToLocal(
               calculationPointOctreeRoot,
               nodesToVisit);
         }
      }
   }
   else
   {
      if(!nodesToVisit.empty())
      {
         for(auto child : _children)
         {
            if(child->isSubdivided() || !child->quadratures().empty())
            {
               auto newNodesToVisit = calculationPointOctreeRoot->getAllNodesAsSet();

               child->translateMultipoleExpansionsToLocal(
                  calculationPointOctreeRoot,
                  newNodesToVisit);
            }
         }
      }
   }
}

std::vector<CalculationPointOctreeNode*> QuadratureOctreeNode::getInteractionList(
   CalculationPointOctreeNode* calculationPointOctreeNode)
{
   std::vector<CalculationPointOctreeNode*> res;

   /*if(calculationPointOctreeNode && 
      calculationPointOctreeNode->parent())
   {
      auto grandparent = calculationPointOctreeNode->parent()->parent();

      if(grandparent)
      {
         auto &parentNearNeighbours = grandparent->children();

         for(auto parentNearNeighbour : parentNearNeighbours)
         {
            for(auto parentNearNeighbourChild : parentNearNeighbour->children())
            {
               if(4 * _box.radius() * _box.radius() < 
                  Vector3::distanceSquared(
                     parentNearNeighbourChild->box().center(),
                     _box.center()))
               {
                  res.push_back(parentNearNeighbourChild);
               }
            }
         }
      }
   }*/

   std::queue< CalculationPointOctreeNode*> queue;
   queue.push(calculationPointOctreeNode);

   while(!queue.empty())
   {
      auto currentNode = queue.front();

      if(4 * _box.radius() * _box.radius() <
         Vector3::distanceSquared(
            currentNode->box().center(),
            _box.center()) && 
         !_box.intersects(currentNode->box()) &&
         !_box.contains(currentNode->box().center()) &&
         !currentNode->box().contains(_box.center()))
      {
         res.push_back(currentNode);
      }
      else
      {
         for(auto child : currentNode->children())
            queue.push(child);
      }
   }

   return res;
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
