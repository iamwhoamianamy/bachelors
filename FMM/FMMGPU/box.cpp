#include <math.h>
#include "box.hpp"

Box::Box()
{
}

Box::Box(const Vector3& center, const Vector3& halfDimensions):
   _center(center), _halfDimensions(halfDimensions)
{
}

bool isInRange(real value, real rangeCenter, real rangeHalfWidth)
{
   return (rangeCenter - rangeHalfWidth <= value &&
           value < rangeCenter + rangeHalfWidth);
}

bool Box::contains(const Vector3& point) const
{
   return (isInRange(point.x, _center.x, _halfDimensions.x) &&
           isInRange(point.y, _center.y, _halfDimensions.y) &&
           isInRange(point.z, _center.z, _halfDimensions.z));
}

bool Box::intersects(const Box& _box) const
{
   return Box(_box.center(), _halfDimensions + _box.halfDimensions()).contains(center());
}

real Box::radius() const
{
    return sqrt(_halfDimensions.x * _halfDimensions.x + 
                _halfDimensions.y * _halfDimensions.y);
}

const Vector3& Box::center() const
{
   return _center;
}

const Vector3& Box::halfDimensions() const
{
   return _halfDimensions;
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
