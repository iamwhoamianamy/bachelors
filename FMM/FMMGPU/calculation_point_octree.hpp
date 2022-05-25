#pragma once
#include <vector>
#include <set>

#include "vector3.cuh"
#include "box.hpp"
#include "harmonic_series.hpp"

class QuadratureOctreeNode;
class CalculationPointOctreeNode;

struct FFMResult
{
   Vector3 point;
   Vector3 answer;
   CalculationPointOctreeNode* node;
};

class CalculationPointOctreeNode
{
private:
   size_t _capacity = 0;
   Box _box;
   std::vector<Vector3*> _points;
   std::vector<CalculationPointOctreeNode*> _children;
   HarmonicSeries<Vector3> _localExpansion;
   CalculationPointOctreeNode* _parent;
 
public:
   //std::set<QuadratureOctreeNode*> leftToInteractWith;

   CalculationPointOctreeNode();
   CalculationPointOctreeNode(
      const Box& box,
      const size_t capacity,
      CalculationPointOctreeNode* parent = nullptr);

   void insert(Vector3& point);
   void insert(std::vector<Vector3>& points);
   void subdivide();

   const Box& box() const;
   bool isSubdivided() const;
   bool isUsefullLeaf() const;
   CalculationPointOctreeNode* parent();
   CalculationPointOctreeNode* parent() const;
   std::vector<CalculationPointOctreeNode*>& children();
   const std::vector<CalculationPointOctreeNode*>& children() const;
   const std::vector<Vector3*>& points() const;
   HarmonicSeries<Vector3>& localExpansion();
   const HarmonicSeries<Vector3>& localExpansion() const;

   std::vector<CalculationPointOctreeNode*> getAllNodes();
   std::set<CalculationPointOctreeNode*> getAllNodesAsSet();
   std::vector<Vector3*> getAllPoints() const;
   size_t getAllNodeCount() const;
   std::vector<FFMResult> calcA(size_t pointCount) ;

   void initAllLocalExpansions(size_t order);
   void propagateLocalExpansions() const;

   ~CalculationPointOctreeNode();
private:
   void calcA(std::vector<FFMResult>& result);

};