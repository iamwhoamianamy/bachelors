#pragma once
#include "real.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

namespace triangle_quadratures
{
   class Vector3
   {
   public:
      real x;
      real y;
      real z;

      __device__ __host__ Vector3(real x = 0, real y = 0, real z = 0);

      __device__ __host__ real& operator[](int i);

      __device__ __host__ Vector3 operator +(const Vector3& rhs) const;
      __device__ __host__ Vector3 operator -(const Vector3& rhs) const;
      __device__ __host__ Vector3 operator *(const real fac) const;
      __device__ __host__ Vector3 operator /(const real fac) const;

      __device__ __host__ Vector3 operator +=(const Vector3& rhs);
      __device__ __host__ Vector3 operator -=(const Vector3& rhs);
      __device__ __host__ Vector3& operator *=(const real fac);
      __device__ __host__ Vector3& operator /=(const real fac);

      __device__ __host__ real LengthSquared() const;
      __device__ __host__ real Length() const;
      __device__ __host__ Vector3 Normalized() const;
      __device__ __host__ Vector3 Perp() const;

      __device__ __host__ void Normalize();
      __device__ __host__ void Limit(const real limitLength);
      __device__ __host__ void SetLength(const real newLength);
      __device__ __host__ Vector3 LeadingCos();

      __device__ __host__ static Vector3 Direction(const Vector3& from, const Vector3& to);
      __device__ __host__ static real DistanceSquared(const Vector3& vec1, const Vector3& vec2);
   };
}