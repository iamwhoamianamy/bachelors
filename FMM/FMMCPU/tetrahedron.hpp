#pragma once
#include "vector3.hpp"
#include <vector>

class Tetrahedron
{
public:
   Tetrahedron();
   Tetrahedron(Vector3& a, Vector3& b, Vector3& c, Vector3& d);
   std::vector<Vector3> points;
   real Volume() const;

   Vector3& a();
   Vector3& b();
   Vector3& c();
   Vector3& d();

   Vector3 a() const;
   Vector3 b() const;
   Vector3 c() const;
   Vector3 d() const;
};