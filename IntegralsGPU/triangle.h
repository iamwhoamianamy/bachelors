#pragma once
#include <vector>

#include "vector3.cuh"
#include"real.h"

namespace triangle_quadratures
{
   class Triangle
   {
   private:
      real area;
      real CalcArea();

   public:
      Vector3 a;
      Vector3 b;
      Vector3 c;

      Triangle();

      Triangle(real ax, real ay, real az,
               real bx, real by, real bz,
               real cx, real cy, real cz);
      
      Triangle(real ax, real ay,
               real bx, real by,
               real cx, real cy);

      Triangle(Vector3 a, Vector3 b, Vector3 c);

      real Area() const;

      std::vector<real> Xs() const;
      std::vector<real> Ys() const;
      std::vector<real> Zs() const;

      Vector3& operator[](int i);

      real XFromST(real s, real t) const;
      real YFromST(real s, real t) const;
      real ZFromST(real s, real t) const;
      Vector3 PointFromST(real s, real t) const;
      Vector3 Normal() const;
   };
}
