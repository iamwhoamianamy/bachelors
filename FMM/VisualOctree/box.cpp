#include "box.hpp"
#include <math.h>

Box::Box(const Vector3& center, const Vector3& halfDimensions):
   center(center), halfDimensions(halfDimensions)
{
}

bool Box::contains(const Vector3& point) const
{
   return (fabs(center.x - point.x) <= halfDimensions.x) &&
          (fabs(center.y - point.y) <= halfDimensions.y) &&
          (fabs(center.z - point.z) <= halfDimensions.z);
}

bool Box::intersects(const Box& _box) const
{
   return Box(_box.center, halfDimensions + _box.halfDimensions).contains(center);
}
//
//real Box::left() const
//{
//   return center.x - halfDimensions.x;
//}
//
//real Box::right() const
//{
//   return center.x + halfDimensions.x;
//}
//
//real Box::top() const
//{
//   return center.y - halfDimensions.y;
//}
//
//real Box::bot() const
//{
//   return center.y + halfDimensions.y;
//}
//
//real Box::near() const
//{
//   return center.z - halfDimensions.z;
//}
//
//real Box::far() const
//{
//   return center.z + halfDimensions.z;
//}
