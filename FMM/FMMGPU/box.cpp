#include <math.h>
#include "box.hpp"

Box::Box()
{
}

Box::Box(const Vector3& center, const Vector3& halfDimensions):
   center(center), halfDimensions(halfDimensions)
{
}

bool isInRange(real value, real rangeCenter, real rangeHalfWidth)
{
   return (rangeCenter - rangeHalfWidth <= value &&
           value < rangeCenter + rangeHalfWidth);
}

bool Box::contains(const Vector3& point) const
{
   return (isInRange(point.x, center.x, halfDimensions.x) &&
           isInRange(point.y, center.y, halfDimensions.y) &&
           isInRange(point.z, center.z, halfDimensions.z));
}

bool Box::intersects(const Box& _box) const
{
   return Box(_box.center, halfDimensions + _box.halfDimensions).contains(center);
}

real Box::radius() const
{
    return sqrt(halfDimensions.x * halfDimensions.x + 
                halfDimensions.y * halfDimensions.y);
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
