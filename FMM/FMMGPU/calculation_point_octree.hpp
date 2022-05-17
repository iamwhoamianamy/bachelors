#pragma once
#include <vector>

#include "vector3.cuh"
#include "box.hpp"
#include "quadrature.hpp"
#include "harmonic_series.hpp"

class CalculationPointOctreeNode
{
private:
   size_t _capacity = 0;
   Box _box;
   std::vector<Vector3*> _points;
   std::vector<CalculationPointOctreeNode*> _children;
   HarmonicSeries<Vector3> _multipoleExpansion;
   HarmonicSeries<Vector3> _localExpansion;
   CalculationPointOctreeNode* _parent;

public:
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
   CalculationPointOctreeNode* parent();
   CalculationPointOctreeNode* parent() const;
   std::vector<CalculationPointOctreeNode*>& children();
   const std::vector<CalculationPointOctreeNode*>& children() const;
   std::vector<Vector3*> points() const;

   std::vector<Vector3*> getAllPoints() const;
   size_t getAllNodeCount() const;

   ~CalculationPointOctreeNode();
};