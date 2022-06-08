#include "reference_cylinder_data.hpp"

ReferenceCylinderData::ReferenceCylinderData(
   size_t id,
   const Vector3& point,
   const Vector3& B,
   real lengthOfB) :
   id(id), point(point), B(B), lengthOfB(lengthOfB)
{
}