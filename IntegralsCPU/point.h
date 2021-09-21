#pragma once
namespace triangle_quadratures
{
   class Point
   {
   public:
      double x;
      double y;
      double z;

      Point(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

      double& operator[](int i);
   };
}