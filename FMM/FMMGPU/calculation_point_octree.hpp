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
   std::vector<Quadrature*> _quadratures;
   std::vector<CalculationPointOctreeNode*> _children;
   HarmonicSeries<Vector3> _multipoleExpansion;
   CalculationPointOctreeNode* _parent;

public:
   CalculationPointOctreeNode();
   CalculationPointOctreeNode(
      const Box& box,
      const size_t capacity,
      CalculationPointOctreeNode* parent = nullptr);

   void insert(Quadrature& point);
   void insert(std::vector<Quadrature>& points);
   void subdivide();
   void quarry(const Box& box, std::vector<Quadrature*>& found);

   const Box& box() const;
   bool isSubdivided() const;
   CalculationPointOctreeNode* parent();
   CalculationPointOctreeNode* parent() const;
   std::vector<CalculationPointOctreeNode*>& children();
   const std::vector<CalculationPointOctreeNode*>& children() const;
   std::vector<Quadrature*> quadratures() const;
   HarmonicSeries<Vector3>& multipoleExpansion();
   const HarmonicSeries<Vector3>& multipoleExpansion() const;

   std::vector<Quadrature*> getAllQuadratures() const;
   void calcLocalMultipolesWithoutTranslation(int n);
   void calcLocalMultipolesWithComplexTranslation(int n);
   void calcLocalMultipolesWithRealTranslation(int n);
   void initAllMultipoleExpansions(size_t n);
   void calcLocalMultipolesAtLeaves(size_t n);
   Vector3 calcA(const Vector3& point) const;
   Vector3 caclRot(const Vector3& point) const;
   size_t getAllNodeCount() const;

   ~CalculationPointOctreeNode();
private:
   bool isLeafAndUseful() const;
};