#pragma once
#include "vector3.hpp"

struct Box
{
   Vector3 center;
   Vector3 halfDimensions;

   Box(const Vector3& center, const Vector3& halfDimensions);
   bool contains(const Vector3& point) const;
   bool intersects(const Box& _box) const;
};