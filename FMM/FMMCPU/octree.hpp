#pragma once
#include <vector>

#include "vector3.hpp"
#include "box.hpp"

class Octree
{
private:
   bool _isSubdivided;
   size_t _capacity;
   Box _box;
   std::vector<Vector3*> _points;
   std::vector<Octree*> _children;

public:
   Octree(const Box& _box, const size_t capacity);

   void insert(Vector3& point);
   void insert(std::vector<Vector3>& points);
   void subdivide();
   void quarry(const Box& box, std::vector<Vector3*>& found);

   const Box& box() const;
   const bool isSubdivided() const;
   const std::vector<Octree*> children() const;

   ~Octree();
};