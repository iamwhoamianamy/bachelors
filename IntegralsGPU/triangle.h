#pragma once
#include <vector>
#include "vector3.cuh"

namespace triangle_quadratures
{
   class Triangle
   {
   private:
      float area;
      float CalcArea();

   public:
      Vector3 a;
      Vector3 b;
      Vector3 c;

      Triangle();

      Triangle(float ax, float ay, float az,
               float bx, float by, float bz,
               float cx, float cy, float cz);
      
      Triangle(float ax, float ay,
               float bx, float by,
               float cx, float cy);

      Triangle(Vector3 a, Vector3 b, Vector3 c);

      float Area() const;

      std::vector<float> Xs() const;
      std::vector<float> Ys() const;
      std::vector<float> Zs() const;

      Vector3& operator[](int i);

      float XFromST(float s, float t) const;
      float YFromST(float s, float t) const;
      float ZFromST(float s, float t) const;
      Vector3 PointFromST(float s, float t) const;
      Vector3 Normal() const;
   };
}
