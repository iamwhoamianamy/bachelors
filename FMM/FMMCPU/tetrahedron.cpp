#include "tetrahedron.hpp"

Tetrahedron::Tetrahedron()
{
   points.resize(4);
}

Tetrahedron::Tetrahedron(Vector3& a, Vector3& b, Vector3& c, Vector3& d)
   :points({a, b, c, d})
{
}

real Tetrahedron::Volume() const
{
   return abs(Vector3::Dot(a() - d(), Vector3::Cross(b() - d(), c() - d()))) / 6.0;
}

Vector3& Tetrahedron::a()
{
   return points[0];
}

Vector3 Tetrahedron::a() const
{
   return points[0];
}

Vector3& Tetrahedron::b()
{
   return points[1];
}

Vector3 Tetrahedron::b() const
{
   return points[1];
}

Vector3& Tetrahedron::c()
{
   return points[2];
}

Vector3 Tetrahedron::c() const
{
   return points[2];
}

Vector3& Tetrahedron::d()
{
   return points[3];
}

Vector3 Tetrahedron::d() const
{
   return points[3];
}
