#pragma once
#include "vector3.cuh"

class Box
{
private:
   Vector3 _center;
   Vector3 _halfDimensions;

public:
   Box();
   Box(const Vector3& center, const Vector3& halfDimensions);
   bool contains(const Vector3& point) const;
   bool contains(const Box& box) const;
   bool intersects(const Box& _box) const;
   real radius() const;

   const Vector3& center() const;
   const Vector3& halfDimensions() const;
};