#pragma once
#include <ostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "real.hpp"

class Vector3
{
public:
   real x;
   real y;
   real z;

   __all__ Vector3(const Vector3& vector);
   __all__ Vector3(Vector3&& vector) noexcept;
   __all__ Vector3(real x = 0, real y = 0, real z = 0);

   __all__ real& operator[](int i);

   __all__ Vector3& operator=(Vector3&& vector) noexcept;
   __all__ Vector3& operator=(const Vector3& vector);

   __all__ Vector3 operator +(const Vector3& rhs) const;
   __all__ Vector3 operator -(const Vector3& rhs) const;
   __all__ Vector3 operator *(const real fac) const;
   __all__ Vector3 operator /(const real fac) const;

   __all__ Vector3& operator +=(const Vector3& rhs);
   __all__ Vector3& operator -=(const Vector3& rhs);
   __all__ Vector3& operator *=(const real fac);
   __all__ Vector3& operator /=(const real fac);

   __all__ real lengthSquared() const;
   __all__ real length() const;
   __all__ Vector3 normalized() const;
   __all__ Vector3 perp() const;

   __all__ void normalize();
   __all__ void limit(const real limitLength);
   __all__ void setLength(const real newLength);
   __all__ Vector3 leadingCos();

   __all__ real r() const;
   __all__ real t() const;
   __all__ real f() const;

   __all__ static Vector3 direction(const Vector3& from, const Vector3& to);
   __all__ static real distanceSquared(const Vector3& vec1, const Vector3& vec2);
   __all__ static real dot(const Vector3& vec1, const Vector3& vec2);
   __all__ static Vector3 cross(const Vector3& vec1, const Vector3& vec2);

   __all__ static Vector3 xAxis();
   __all__ static Vector3 yAxis();
   __all__ static Vector3 zAxis();

   __all__ void printWithWidth(std::ostream& os, size_t width);
};

std::ostream& operator<<(std::ostream& os, const Vector3& vec);