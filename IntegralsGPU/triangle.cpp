#include <cmath>
#include "triangle.h"
#include "exeptions.h"

namespace triangle_quadratures
{
   real Triangle::CalcArea()
   {
      //return 0.5 * abs((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
      Vector3 ab(b.x - a.x, b.y - a.y, b.z - a.z);
      Vector3 ac(c.x - a.x, c.y - a.y, c.z - a.z);

      real t0 = ab.y * ac.z - ab.z * ac.y;
      real t1 = ab.z * ac.x - ab.x * ac.z;
      real t2 = ab.x * ac.y - ab.y * ac.x;

      return 0.5 * sqrt(t0 * t0 + t1 * t1 + t2 * t2);
   }

   Triangle::Triangle() : a(), b(), c()
   {
      area = 0;
   }

   Triangle::Triangle(real ax, real ay, real az,
                      real bx, real by, real bz,
                      real cx, real cy, real cz) : 
      a(ax, ay, az), b(bx, by, bz), c(cx, cy, cz)
   {
      area = CalcArea();
   }

   Triangle::Triangle(real ax, real ay,
                      real bx, real by,
                      real cx, real cy) :
      a(ax, ay, 0), b(bx, by, 0), c(cx, cy, 0)
   {
      area = CalcArea();
   }

   Triangle::Triangle(Vector3 a, Vector3 b, Vector3 c) : a(a), b(b), c(c)
   {
      area = CalcArea();
   }

   real Triangle::Area() const { return area; }

   std::vector<real> Triangle::Xs() const
   {
      return { a.x, b.x, c.x };
   }

   std::vector<real> Triangle::Ys() const
   {
      return { a.y, b.y, c.y };
   }

   std::vector<real> Triangle::Zs() const
   {
      return { a.z, b.z, c.z };
   }

   Vector3& Triangle::operator[](int i)
   {
      switch(i)
      {
         case 0: return a;
         case 1: return b;
         case 2: return c;
         default: throw RangeExeption();
      }
   }

   real Triangle::XFromST(real s, real t) const
   {
      //if(s > 1.0 || t > 1.0 || t > 1.0 - s)
      //   throw RangeExeption();

      return a.x + (b.x - a.x) * s + (c.x - a.x) * t;
   }

   real Triangle::YFromST(real s, real t) const
   {
      //if(s > 1.0 || t > 1.0 || t > 1.0 - s)
      //   throw RangeExeption();

      return a.y + (b.y - a.y) * s + (c.y - a.y) * t;
   }

   real Triangle::ZFromST(real s, real t) const
   {
      //if(s > 1.0 || t > 1.0 || t > 1.0 - s)
      //   throw RangeExeption();

      return a.z + (b.z - a.z) * s + (c.z - a.z) * t;
   }

   Vector3 Triangle::PointFromST(real s, real t) const
   {
      //if(s > 1.0 || t > 1.0 || t > 1.0 - s)
      //   throw RangeExeption();

      return Vector3(XFromST(s, t), YFromST(s, t), ZFromST(s, t));
   }

   Vector3 Triangle::Normal() const
   {
      Vector3 u = a - b;
      Vector3 v = a - c;

      Vector3 normal(
         u.y * v.z - u.z * v.y,
         u.z * v.x - u.x * v.z,
         u.x * v.y - u.y * v.x
      );

      return normal.Normalized();
   }
}