#pragma once
#include <vector>

#include "vector3.cuh"

class Tetrahedron
{
public:
   Tetrahedron();
   Tetrahedron(Vector3& a, Vector3& b, Vector3& c, Vector3& d);
   std::vector<Vector3> points;
   real volume() const;

   Vector3& a();
   Vector3& b();
   Vector3& c();
   Vector3& d();

   Vector3 a() const;
   Vector3 b() const;
   Vector3 c() const;
   Vector3 d() const;
};