#pragma once
namespace triangle_quadratures
{
   class Vector3
   {
   public:
      double x;
      double y;
      double z;

      Vector3(double x = 0, double y = 0, double z = 0);

      double& operator[](int i);

      Vector3 operator +(const Vector3& rhs) const;
      Vector3 operator -(const Vector3& rhs) const;
      Vector3 operator *(const double fac) const;
      Vector3 operator /(const double fac) const;

      Vector3 operator +=(const Vector3& rhs);
      Vector3 operator -=(const Vector3& rhs);
      Vector3& operator *=(const double fac);
      Vector3& operator /=(const double fac);

      double LengthSquared() const;
      double Length() const;
      Vector3 Normalized() const;
      Vector3 Perp() const;

      void Normalize();
      void Limit(const double limitLength);
      void SetLength(const double newLength);
      Vector3 LeadingCos();

      static Vector3 Direction(const Vector3& from, const Vector3& to);
      static double DistanceSquared(const Vector3& vec1, const Vector3& vec2);
   };
}