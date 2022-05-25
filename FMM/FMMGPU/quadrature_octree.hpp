#pragma once
#include <vector>
#include <set>

#include "vector3.cuh"
#include "box.hpp"
#include "quadrature.hpp"
#include "harmonic_series.hpp"
#include "calculation_point_octree.hpp"

class QuadratureOctreeNode
{
private:
   size_t _capacity = 0;
   Box _box;
   std::vector<Quadrature*> _quadratures;
   std::vector<QuadratureOctreeNode*> _children;
   HarmonicSeries<Vector3> _multipoleExpansion;
   QuadratureOctreeNode* _parent;

public:

   QuadratureOctreeNode();
   QuadratureOctreeNode(
      const Box& box,
      const size_t capacity,
      QuadratureOctreeNode* parent = nullptr);

   void insert(Quadrature& point);
   void insert(std::vector<Quadrature>& points);
   void subdivide();

   const Box& box() const;
   bool isSubdivided() const;
   bool isUsefullLeaf() const;
   QuadratureOctreeNode* parent();
   QuadratureOctreeNode* parent() const;
   std::vector<QuadratureOctreeNode*>& children();
   const std::vector<QuadratureOctreeNode*>& children() const;
   const std::vector<Quadrature*>& quadratures() const;
   HarmonicSeries<Vector3>& multipoleExpansion();
   const HarmonicSeries<Vector3>& multipoleExpansion() const;

   std::vector<Quadrature*> getAllQuadratures() const;

   void calcMultipoleExpansionsWithoutTranslation(int n);
   void calcMultipoleExpansionsWithComplexTranslation(int n);
   void calcMultipoleExpansionsWithRealTranslation(int n);
   void calcMultipoleExpansionsAtLeaves(size_t n);
   void initAllMultipoleExpansions(size_t n);

   Vector3 calcA(const Vector3& point) const;
   Vector3 caclRot(const Vector3& point) const;
   size_t getAllNodeCount() const;

   ~QuadratureOctreeNode();

private:
   bool isLeafAndUseful() const;
};