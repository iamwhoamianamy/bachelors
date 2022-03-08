#pragma once
#include "../FMMCPU/vector3.hpp"

struct Box
{
   Vector3 center;
   Vector3 halfDimensions;

   Box(const Vector3& center, const Vector3& halfDimensions);
   bool doContain(const Vector3& point) const;
   bool doIntersect(const Box& _box) const;
};