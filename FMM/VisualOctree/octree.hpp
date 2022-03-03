#pragma once
#include <vector>

#include "../FMMCPU/vector3.hpp"
#include "box.hpp"

class Octree
{
private:
   std::vector<Octree*> _children;
   Box _box;
   bool _isSubdivided;
   size_t _capacity;
   std::vector<Vector3*> points;

public:
   Octree(const Box& _box, const size_t capacity);
   void insert(Vector3& point);
   void insert(std::vector<Vector3>& points);
   void subdivide();
   void quarry(const Box& box, std::vector<Vector3*>& found);

   const Box& box() const;
   const bool isSubdivided() const;
   const std::vector<Octree*> children() const;

   /*Octree* botSouthWest();
   Octree* botSouthEast();
   Octree* botNorthEast();
   Octree* botNorthWest();

   Octree* topSouthWest();
   Octree* topSouthEast();
   Octree* topNorthEast();
   Octree* topNorthWest();*/

};