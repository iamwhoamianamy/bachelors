#pragma once
#include "vector3.cuh"

struct ReferenceCylinderData
{
   size_t id;
   Vector3 point;
   Vector3 B;
   real lengthOfB;

   ReferenceCylinderData(
      size_t id,
      const Vector3& point,
      const Vector3& B,
      real lengthOfB);
};
