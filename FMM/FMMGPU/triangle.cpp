#include "triangle.hpp"

Triangle::Triangle()
{
   points.resize(3);
}

Triangle::Triangle(const Vector3& a, const Vector3& b, const Vector3& c)
   :points({a, b, c})
{
}

real Triangle::square() const
{
   Vector3 ab(b().x - a().x, b().y - a().y, b().z - a().z);
   Vector3 ac(c().x - a().x, c().y - a().y, c().z - a().z);

   real t0 = ab.y * ac.z - ab.z * ac.y;
   real t1 = ab.z * ac.x - ab.x * ac.z;
   real t2 = ab.x * ac.y - ab.y * ac.x;

   return 0.5 * sqrt(t0 * t0 + t1 * t1 + t2 * t2);
}

Vector3& Triangle::a()
{
   return points[0];
}

Vector3 Triangle::a() const
{
   return points[0];
}

Vector3& Triangle::b()
{
   return points[1];
}

Vector3 Triangle::b() const
{
   return points[1];
}

Vector3& Triangle::c()
{
   return points[2];
}

Vector3 Triangle::c() const
{
   return points[2];
}