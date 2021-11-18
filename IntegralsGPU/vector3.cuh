#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

namespace triangle_quadratures
{
   class Vector3
   {
   public:
      float x;
      float y;
      float z;

      __device__ __host__ Vector3(float x = 0, float y = 0, float z = 0);

      __device__ __host__ float& operator[](int i);

      __device__ __host__ Vector3 operator +(const Vector3& rhs) const;
      __device__ __host__ Vector3 operator -(const Vector3& rhs) const;
      __device__ __host__ Vector3 operator *(const float fac) const;
      __device__ __host__ Vector3 operator /(const float fac) const;

      __device__ __host__ Vector3 operator +=(const Vector3& rhs);
      __device__ __host__ Vector3 operator -=(const Vector3& rhs);
      __device__ __host__ Vector3& operator *=(const float fac);
      __device__ __host__ Vector3& operator /=(const float fac);

      __device__ __host__ float LengthSquared() const;
      __device__ __host__ float Length() const;
      __device__ __host__ Vector3 Normalized() const;
      __device__ __host__ Vector3 Perp() const;

      __device__ __host__ void Normalize();
      __device__ __host__ void Limit(const float limitLength);
      __device__ __host__ void SetLength(const float newLength);
      __device__ __host__ Vector3 LeadingCos();

      __device__ __host__ static Vector3 Direction(const Vector3& from, const Vector3& to);
      __device__ __host__ static float DistanceSquared(const Vector3& vec1, const Vector3& vec2);
   };
}