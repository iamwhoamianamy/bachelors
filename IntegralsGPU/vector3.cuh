#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

namespace triangle_quadratures
{
   class Vector3
   {
   public:
      double x;
      double y;
      double z;

      __device__ __host__ Vector3(double x = 0, double y = 0, double z = 0);

      __device__ __host__ double& operator[](int i);

      __device__ __host__ Vector3 operator +(const Vector3& rhs) const;
      __device__ __host__ Vector3 operator -(const Vector3& rhs) const;
      __device__ __host__ Vector3 operator *(const double fac) const;
      __device__ __host__ Vector3 operator /(const double fac) const;

      __device__ __host__ Vector3 operator +=(const Vector3& rhs);
      __device__ __host__ Vector3 operator -=(const Vector3& rhs);
      __device__ __host__ Vector3& operator *=(const double fac);
      __device__ __host__ Vector3& operator /=(const double fac);

      __device__ __host__ double LengthSquared() const;
      __device__ __host__ double Length() const;
      __device__ __host__ Vector3 Normalized() const;
      __device__ __host__ Vector3 Perp() const;

      __device__ __host__ void Normalize();
      __device__ __host__ void Limit(const double limitLength);
      __device__ __host__ void SetLength(const double newLength);
      __device__ __host__ Vector3 LeadingCos();

      __device__ __host__ static Vector3 Direction(const Vector3& from, const Vector3& to);
      __device__ __host__ static double DistanceSquared(const Vector3& vec1, const Vector3& vec2);
   };
}