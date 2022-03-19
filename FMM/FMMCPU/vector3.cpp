#include "vector3.hpp"
#include <exception>

Vector3::Vector3(const Vector3& vector)
{
   x = vector.x;
   y = vector.y;
   z = vector.z;
}

Vector3::Vector3(Vector3&& vector) noexcept
{
   x = std::move(vector.x);
   y = std::move(vector.y);
   z = std::move(vector.z);
}

Vector3::Vector3(real x, real y, real z)
   : x(x), y(y), z(z) {}

real& Vector3::operator[](int i)
{
   switch(i)
   {
      case 0: return x;
      case 1: return y;
      case 2: return z;
      default: throw std::exception();
   }
}

Vector3& Vector3::operator=(Vector3&& vector) noexcept
{
   x = std::move(vector.x);
   y = std::move(vector.y);
   z = std::move(vector.z);

   return *this;
}

Vector3& Vector3::operator=(const Vector3& vector)
{
   x = vector.x;
   y = vector.y;
   z = vector.z;

   return *this;
}

Vector3 Vector3::operator +(const Vector3& rhs) const
{
   Vector3 res;
   res.x = this->x + rhs.x;
   res.y = this->y + rhs.y;
   res.z = this->z + rhs.z;

   return res;
}

Vector3 Vector3::operator -(const Vector3& rhs) const
{
   Vector3 res;
   res.x = this->x - rhs.x;
   res.y = this->y - rhs.y;
   res.z = this->z - rhs.z;

   return res;
}

Vector3 Vector3::operator *(const real fac) const
{
   Vector3 res;
   res.x = this->x * fac;
   res.y = this->y * fac;
   res.z = this->z * fac;

   return res;
}

Vector3 Vector3::operator /(const real fac) const
{
   Vector3 res;
   res.x = this->x / fac;
   res.y = this->y / fac;
   res.z = this->z / fac;

   return res;
}

Vector3& Vector3::operator +=(const Vector3& rhs)
{
   this->x += rhs.x;
   this->y += rhs.y;
   this->z += rhs.z;

   return *this;
}

Vector3& Vector3::operator -=(const Vector3& rhs)
{
   this->x -= rhs.x;
   this->y -= rhs.y;
   this->z -= rhs.z;

   return *this;
}

Vector3& Vector3::operator *=(const real fac)
{
   this->x *= fac;
   this->y *= fac;
   this->z *= fac;

   return *this;
}

Vector3& Vector3::operator /=(const real fac)
{
   if(fac != 0)
   {
      this->x /= fac;
      this->y /= fac;
      this->z /= fac;
   }
   else
   {
      this->x = 0;
      this->y = 0;
      this->z = 0;
   }

   return *this;
}

real Vector3::lengthSquared() const
{
   return this->x * this->x +
      this->y * this->y +
      this->z * this->z;
}

real Vector3::length() const
{
   return sqrt(this->x * this->x +
               this->y * this->y +
               this->z * this->z);
}

Vector3 Vector3::normalized() const
{
   Vector3 res;
   real length = length();

   if(length)
   {
      res.x = this->x / length;
      res.y = this->y / length;
      res.z = this->z / length;
   }

   return res;
}

void Vector3::normalize()
{
   real length = length();

   if(length != 0)
   {
      this->x /= length;
      this->y /= length;
      this->z /= length;
   }
}

void Vector3::limit(const real limitLength)
{
   real length = length();

   if(length != 0 && length > limitLength)
   {
      this->x = this->x / length * limitLength;
      this->y = this->y / length * limitLength;
      this->z = this->z / length * limitLength;
   }
}

void Vector3::setLength(const real newLength)
{
   normalize();
   this->x *= newLength;
   this->y *= newLength;
   this->z *= newLength;
}

Vector3 Vector3::perp() const
{
   return Vector3(-y, x);
}

Vector3 Vector3::direction(const Vector3& from, const Vector3& to)
{
   Vector3 res = to - from;
   return res.normalized();
}

real Vector3::distanceSquared(const Vector3& vec1, const Vector3& vec2)
{
   return (vec1 - vec2).lengthSquared();
}

double Vector3::dot(const Vector3& vec1, const Vector3& vec2)
{
   return vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
}

Vector3 Vector3::cross(const Vector3& vec1, const Vector3& vec2)
{
   return Vector3(vec1.y * vec2.z - vec1.z * vec2.y,
                  vec1.z * vec2.x - vec1.x * vec2.z,
                  vec1.x * vec2.y - vec1.y * vec2.x);
}

Vector3 Vector3::leadingCos()
{
   real l = length();

   return Vector3(x / l, y / l, z / l);
}

real Vector3::r() const
{
   return sqrt(x * x + y * y + z * z);
}

real Vector3::t() const
{
   return atan(sqrt(x * x + y + y) / z);
}

real Vector3::f() const
{
   return atan2(y, x);
}
