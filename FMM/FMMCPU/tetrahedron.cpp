#include "tetrahedron.hpp"

Tetrahedron::Tetrahedron()
{
   points.resize(4);
}

Tetrahedron::Tetrahedron(Vector3& a, Vector3& b, Vector3& c, Vector3& d)
   :points({a, b, c, d})
{
}

Vector3& Tetrahedron::a()
{
   return points[0];
}

Vector3& Tetrahedron::b()
{
   return points[1];
}

Vector3& Tetrahedron::c()
{
   return points[2];
}

Vector3& Tetrahedron::d()
{
   return points[3];
}
