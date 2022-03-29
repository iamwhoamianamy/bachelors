#pragma once
#include "real.hpp"
#include <ostream>

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

   real lengthSquared() const;
   real length() const;
   Vector3 normalized() const;
   Vector3 perp() const;

   void normalize();
   void limit(const real limitLength);
   void setLength(const real newLength);
   Vector3 leadingCos();

   real r() const;
   real t() const;
   real f() const;

   static Vector3 direction(const Vector3& from, const Vector3& to);
   static real distanceSquared(const Vector3& vec1, const Vector3& vec2);
   static double dot(const Vector3& vec1, const Vector3& vec2);
   static Vector3 cross(const Vector3& vec1, const Vector3& vec2);

   static Vector3 xAxis();
   static Vector3 yAxis();
   static Vector3 zAxis();

   void printWithWidth(std::ostream& os, size_t width);
};

std::ostream& operator<<(std::ostream& os, const Vector3& vec);
