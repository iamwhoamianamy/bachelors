#pragma once
#include <vector>
#include "vector3.h"

namespace triangle_quadratures
{
   class Triangle
   {
   private:
      double area;
      double CalcArea();

   public:
      Vector3 a;
      Vector3 b;
      Vector3 c;

      Triangle();

      Triangle(double ax, double ay, double az,
               double bx, double by, double bz,
               double cx, double cy, double cz);
      
      Triangle(double ax, double ay,
               double bx, double by,
               double cx, double cy);

      Triangle(Vector3 a, Vector3 b, Vector3 c);

      double Area() const;

      std::vector<double> Xs() const;
      std::vector<double> Ys() const;
      std::vector<double> Zs() const;

      Vector3& operator[](int i);

      double XFromST(double s, double t) const;
      double YFromST(double s, double t) const;
      double ZFromST(double s, double t) const;
      Vector3 PointFromST(double s, double t) const;
      Vector3 Normal() const;
   };
}
