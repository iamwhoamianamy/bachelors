#pragma once
#include <vector>

#include "vector3.hpp"
#include "box.hpp"
#include "quadrature.hpp"
#include "harmonic_series.hpp"

class Octree
{
private:
   bool _isSubdivided = false;
   size_t _capacity = 0;
   Box _box;
   std::vector<Quadrature*> _quadratures;
   std::vector<Octree*> _children;
   HarmonicSeries<Vector3> _multipoleExpansion;

public:

   Octree();
   Octree(const Box& box, const size_t capacity);

   void insert(Quadrature& point);
   void insert(std::vector<Quadrature>& points);
   void subdivide();
   void quarry(const Box& box, std::vector<Quadrature*>& found);

   const Box& box() const;
   const bool isSubdivided() const;
   const std::vector<Octree*> children() const;
   const HarmonicSeries<Vector3>& multipoleExpansion() const;

   std::vector<Quadrature*> getAllQuadratures() const;
   void calcLocalMultipolesWithoutTranslation(int n);
   void calcLocalMultipolesWithTranslation(int n);
   Vector3 calcA(const Vector3& point) const;
   Vector3 caclRot(const Vector3& point) const;

   ~Octree();
   
private:
   void accountChildContribution(Octree* child, int n);
};