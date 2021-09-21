#pragma once
#include <vector>
#include "point.h"

namespace triangle_quadratures
{
   struct Triangle
   {
   private:
      double area;
      double CalcArea();

   public:
      Point a;
      Point b;
      Point c;

      Triangle();

      Triangle(double ax, double ay, double az,
               double bx, double by, double bz,
               double cx, double cy, double cz);
      
      Triangle(double ax, double ay,
               double bx, double by,
               double cx, double cy);

      Triangle(Point a, Point b, Point c);

      double Area() const;

      std::vector<double> Xs() const;
      std::vector<double> Ys() const;
      std::vector<double> Zs() const;

      Point& operator[](int i);

      double XFromST(double s, double t) const;
      double YFromST(double s, double t) const;
      double ZFromST(double s, double t) const;
      Point PointFromST(double s, double t) const;
   };
}
