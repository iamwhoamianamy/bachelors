#pragma once
#include <vector>

#include "vector3.cuh"

class Triangle
{
public:
   Triangle();
   Triangle(const Vector3& a, const Vector3& b, const Vector3& c);
   std::vector<Vector3> points;
   real square() const;

   Vector3& a();
   Vector3& b();
   Vector3& c();

   Vector3 a() const;
   Vector3 b() const;
   Vector3 c() const;
};