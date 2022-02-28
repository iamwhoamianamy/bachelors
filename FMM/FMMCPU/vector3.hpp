#pragma once
#include "real.hpp"

class Vector3
{
public:
   real x;
   real y;
   real z;

   Vector3(const Vector3& vector);
   Vector3(Vector3&& vector) noexcept;
   Vector3(real x = 0, real y = 0, real z = 0);

   real& operator[](int i);

   Vector3& operator=(Vector3&& vector) noexcept;
   Vector3& operator=(const Vector3& vector);

   Vector3 operator +(const Vector3& rhs) const;
   Vector3 operator -(const Vector3& rhs) const;
   Vector3 operator *(const real fac) const;
   Vector3 operator /(const real fac) const;

   Vector3& operator +=(const Vector3& rhs);
   Vector3& operator -=(const Vector3& rhs);
   Vector3& operator *=(const real fac);
   Vector3& operator /=(const real fac);

   real LengthSquared() const;
   real Length() const;
   Vector3 Normalized() const;
   Vector3 Perp() const;

   void Normalize();
   void Limit(const real limitLength);
   void SetLength(const real newLength);
   Vector3 LeadingCos();

   static Vector3 Direction(const Vector3& from, const Vector3& to);
   static real DistanceSquared(const Vector3& vec1, const Vector3& vec2);
   static double Dot(const Vector3& vec1, const Vector3& vec2);
   static Vector3 Cross(const Vector3& vec1, const Vector3& vec2);
};