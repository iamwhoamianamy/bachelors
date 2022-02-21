#pragma once
#include "vector3.hpp"
#include <vector>

class Tetrahedron
{
private:
   std::vector<Vector3> points;
public:
   Tetrahedron();
   Tetrahedron(Vector3& a, Vector3& b, Vector3& c, Vector3& d);

   Vector3& a();
   Vector3& b();
   Vector3& c();
   Vector3& d();
};