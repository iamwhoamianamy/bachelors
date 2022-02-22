#pragma once
#include "vector3.hpp"
#include <vector>

class Tetrahedron
{
public:
   Tetrahedron();
   Tetrahedron(Vector3& a, Vector3& b, Vector3& c, Vector3& d);
   std::vector<Vector3> points;

   Vector3& a();
   Vector3& b();
   Vector3& c();
   Vector3& d();
};