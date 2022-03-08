#include "octree.hpp"

Octree::Octree(const Box& _box, const size_t capacity) :
   _box(_box), _capacity(capacity), _isSubdivided(false)
{

}

void Octree::insert(std::vector<Vector3>& points)
{
   for(auto& point : points)
   {
      insert(point);
   }
}

void Octree::insert(Vector3& point)
{
   if(!_box.doContain(point))
   {
      return;
   }

   if(points.size() < _capacity)
   {
      if(!_isSubdivided)
      {
         points.push_back(&point);
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
      _isSubdivided = true;

      points.push_back(&point);

      for(auto p : points)
      {
         for(auto& child : _children)
         {
            child->insert(*p);
         }
      }

      points.clear();
   }
}

void Octree::subdivide()
{
   float x = _box.center.x;
   float y = _box.center.y;
   float z = _box.center.z;

   float w = _box.halfDimensions.x;
   float h = _box.halfDimensions.y;
   float d = _box.halfDimensions.z;

   Vector3 childrenHalfDimensions = { w / 2, h / 2, d / 2 };

   _children.reserve(8);

   _children.emplace_back(new Octree(Box({ x - w / 2, y - h / 2, z - d / 2 }, childrenHalfDimensions), _capacity));
   _children.emplace_back(new Octree(Box({ x + w / 2, y - h / 2, z - d / 2 }, childrenHalfDimensions), _capacity));
   _children.emplace_back(new Octree(Box({ x + w / 2, y + h / 2, z - d / 2 }, childrenHalfDimensions), _capacity));
   _children.emplace_back(new Octree(Box({ x - w / 2, y + h / 2, z - d / 2 }, childrenHalfDimensions), _capacity));

   _children.emplace_back(new Octree(Box({ x - w / 2, y - h / 2, z + d / 2 }, childrenHalfDimensions), _capacity));
   _children.emplace_back(new Octree(Box({ x + w / 2, y - h / 2, z + d / 2 }, childrenHalfDimensions), _capacity));
   _children.emplace_back(new Octree(Box({ x + w / 2, y + h / 2, z + d / 2 }, childrenHalfDimensions), _capacity));
   _children.emplace_back(new Octree(Box({ x - w / 2, y + h / 2, z + d / 2 }, childrenHalfDimensions), _capacity));
}

void Octree::quarry(const Box& range, std::vector<Vector3*>& found)
{
   if(!this->_box.doIntersect(range))
   {
      return;
   }
   else
   {
      for(auto point : points)
      {
         if(range.doContain(*point))
         {
            found.push_back(point);
         }
      }

      if(_isSubdivided)
      {
         for(auto& child : _children)
         {
            child->quarry(range, found);
         }
      }
   }
}

const Box& Octree::box() const
{
   return this->_box;
}

const bool Octree::isSubdivided() const
{
   return _isSubdivided;
}

const std::vector<Octree*> Octree::children() const
{
   return _children;
}
