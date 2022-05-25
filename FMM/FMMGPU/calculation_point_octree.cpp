#include "calculation_point_octree.hpp"
#include "quadrature_octree.hpp"

#include <iostream>

#include "integration.hpp"
#include "harmonics.hpp"
#include "multipole_translator.hpp"

#include <thrust/complex.h>

CalculationPointOctreeNode::CalculationPointOctreeNode() : 
   _parent(nullptr)
{
}

CalculationPointOctreeNode::CalculationPointOctreeNode(
   const Box& box,
   const size_t capacity,
   CalculationPointOctreeNode* parent) :
   _parent(parent), _box(box), _capacity(capacity)
{

}

void CalculationPointOctreeNode::insert(std::vector<Vector3>& points)
{
   for(auto& point : points)
   {
      insert(point);
   }
}

void CalculationPointOctreeNode::insert(Vector3& point)
{
   if(!_box.contains(point))
   {
      return;
   }

   if(_points.size() < _capacity)
   {
      if(!isSubdivided())
      {
         _points.push_back(&point);
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

      _points.push_back(&point);

      for(auto p : _points)
      {
         for(auto& child : _children)
         {
            child->insert(*p);
         }
      }

      _points.resize(0);
   }
}

void CalculationPointOctreeNode::subdivide()
{
   float x = _box.center().x;
   float y = _box.center().y;
   float z = _box.center().z;

   float w = _box.halfDimensions().x;
   float h = _box.halfDimensions().y;
   float d = _box.halfDimensions().z;

   Vector3 childrenHalfDimensions = { w / 2, h / 2, d / 2 };

   _children.reserve(8);

   _children.emplace_back(new CalculationPointOctreeNode(Box({ x - w / 2, y - h / 2, z - d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new CalculationPointOctreeNode(Box({ x + w / 2, y - h / 2, z - d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new CalculationPointOctreeNode(Box({ x + w / 2, y + h / 2, z - d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new CalculationPointOctreeNode(Box({ x - w / 2, y + h / 2, z - d / 2 }, childrenHalfDimensions), _capacity, this));

   _children.emplace_back(new CalculationPointOctreeNode(Box({ x - w / 2, y - h / 2, z + d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new CalculationPointOctreeNode(Box({ x + w / 2, y - h / 2, z + d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new CalculationPointOctreeNode(Box({ x + w / 2, y + h / 2, z + d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new CalculationPointOctreeNode(Box({ x - w / 2, y + h / 2, z + d / 2 }, childrenHalfDimensions), _capacity, this));
}

std::vector<CalculationPointOctreeNode*> CalculationPointOctreeNode::getAllNodes()
{
   std::vector<CalculationPointOctreeNode*> result;
   result.push_back(this);

   for(auto child : _children)
   {
      auto childChildren = child->getAllNodes();
      result.insert(result.end(), childChildren.begin(), childChildren.end());
   }

   return result;
}

std::set<CalculationPointOctreeNode*> CalculationPointOctreeNode::getAllNodesAsSet()
{
   std::set<CalculationPointOctreeNode*> result;
   result.insert(this);

   for(auto child : _children)
   {
      if(child->isSubdivided() || !child->points().empty())
      {
         auto childChildren = child->getAllNodesAsSet();
         result.insert(childChildren.begin(), childChildren.end());
      }
   }

   return result;
}

std::vector<Vector3*> CalculationPointOctreeNode::getAllPoints() const
{
   std::vector<Vector3*> result;
   result.insert(result.end(), _points.begin(), _points.end());

   for(auto child : _children)
   {
      auto childPoints = child->getAllPoints();
      result.insert(result.end(), childPoints.begin(), childPoints.end());
   }

   return result;
}

size_t CalculationPointOctreeNode::getAllNodeCount() const
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

std::vector<FFMResult> CalculationPointOctreeNode::calcA(size_t pointCount)
{
   std::vector<FFMResult> res;
   res.reserve(pointCount);
   calcA(res);

   return res;
}

void CalculationPointOctreeNode::initAllLocalExpansions(size_t order)
{
   if(isSubdivided() || !_points.empty())
   {
      _localExpansion = HarmonicSeries<Vector3>(order);

      for(auto child : _children)
      {
         child->initAllLocalExpansions(order);
      }
   }
}

void CalculationPointOctreeNode::calcA(std::vector<FFMResult>& result)
{
   if(!_points.empty())
   {
      for(auto point : _points)
      {
         auto translation = *point - _box.center();
         auto regularHarmonic = Harmonics::calcRegularSolidHarmonics(
            _localExpansion.order(), translation);
         Vector3 aInPoint = mult(_localExpansion, regularHarmonic);
         result.emplace_back(*point, aInPoint , this);
      }
   }
   else
   {
      for (auto child : _children)
      {
         child->calcA(result);
      }
   }
}

void CalculationPointOctreeNode::propagateLocalExpansions() const
{
   if(_localExpansion.order() != 0)
   {
      for(auto child : _children)
      {
         if(child->localExpansion().order() !=0 &&
            (child->isSubdivided() || !child->points().empty()))
         {
            auto translation = _box.center() - child->box().center();
            child->_localExpansion.add(
               MultipoleTranslator::translateLocalWithComplex(
                  _localExpansion,
                  translation));
         }
      }
   }

   for(auto child : _children)
   {
      child->propagateLocalExpansions();
   }
}

const Box& CalculationPointOctreeNode::box() const
{
   return this->_box;
}

bool CalculationPointOctreeNode::isSubdivided() const
{
   return !_children.empty();
}

bool CalculationPointOctreeNode::isUsefullLeaf() const
{
   return !_points.empty();
}

CalculationPointOctreeNode* CalculationPointOctreeNode::parent()
{
   return _parent;
}

CalculationPointOctreeNode* CalculationPointOctreeNode::parent() const
{
   return _parent;
}

std::vector<CalculationPointOctreeNode*>& CalculationPointOctreeNode::children()
{
   return _children;
}

const std::vector<CalculationPointOctreeNode*>& CalculationPointOctreeNode::children() const
{
   return _children;
}

const std::vector<Vector3*>& CalculationPointOctreeNode::points() const
{
   return _points;
}

HarmonicSeries<Vector3>& CalculationPointOctreeNode::localExpansion()
{
   return _localExpansion;
}

const HarmonicSeries<Vector3>& CalculationPointOctreeNode::localExpansion() const
{
   return _localExpansion;
}

CalculationPointOctreeNode::~CalculationPointOctreeNode()
{
   for(auto child : _children)
   {
      delete child;
   }
}