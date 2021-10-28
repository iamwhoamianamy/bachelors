#include "vector3.h"
#include "exeptions.h"

namespace triangle_quadratures
{
   Vector3::Vector3(double x, double y, double z)
      : x(x), y(y), z(z) {}

   double& Vector3::operator[](int i)
   {
      switch(i)
      {
         case 0: return x;
         case 1: return y;
         case 2: return z;
         default: throw RangeExeption();
      }
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

   Vector3 Vector3::operator *(const double fac) const
   {
      Vector3 res;
      res.x = this->x * fac;
      res.y = this->y * fac;
      res.z = this->z * fac;

      return res;
   }

   Vector3 Vector3::operator /(const double fac) const
   {
      Vector3 res;
      res.x = this->x / fac;
      res.y = this->y / fac;
      res.z = this->z / fac;

      return res;
   }

   Vector3 Vector3::operator +=(const Vector3& rhs)
   {
      this->x += rhs.x;
      this->y += rhs.y;
      this->z += rhs.z;

      return *this;
   }

   Vector3 Vector3::operator -=(const Vector3& rhs)
   {
      this->x -= rhs.x;
      this->y -= rhs.y;
      this->z -= rhs.z;

      return *this;
   }

   Vector3& Vector3::operator *=(const double fac)
   {
      this->x *= fac;
      this->y *= fac;
      this->z *= fac;

      return *this;
   }

   Vector3& Vector3::operator /=(const double fac)
   {
      if (fac != 0)
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

   double Vector3::LengthSquared() const
   {
      return this->x * this->x +
             this->y * this->y +
             this->z * this->z;
   }

   double Vector3::Length() const
   {
      return sqrt(this->x * this->x + 
                  this->y * this->y +
                  this->z * this->z);
   }

   Vector3 Vector3::Normalized() const
   {
      Vector3 res;
      double length = Length();

      if (length)
      {
         res.x = this->x / length;
         res.y = this->y / length;
         res.z = this->z / length;
      }

      return res;
   }

   void Vector3::Normalize()
   {
      double length = Length();

      if (length != 0)
      {
         this->x /= length;
         this->y /= length;
         this->z /= length;
      }
   }

   void Vector3::Limit(const double limitLength)
   {
      double length = Length();

      if (length != 0 && length > limitLength)
      {
         this->x = this->x / length * limitLength;
         this->y = this->y / length * limitLength;
         this->z = this->z / length * limitLength;
      }
   }

   void Vector3::SetLength(const double newLength)
   {
      Normalize();
      this->x *= newLength;
      this->y *= newLength;
      this->z *= newLength;
   }

   Vector3 Vector3::Perp() const
   {
      return Vector3(-y, x);
   }

   Vector3 Vector3::Direction(const Vector3& from, const Vector3& to)
   {
      Vector3 res = to - from;
      return res.Normalized();
   }

   double Vector3::DistanceSquared(const Vector3& vec1, const Vector3& vec2)
   {
      return (vec1 - vec2).LengthSquared();
   }

   Vector3 Vector3::LeadingCos()
   {
      double l = Length();

      return Vector3(x / l, y / l, z / l);
   }
}