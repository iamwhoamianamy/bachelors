#include <iostream>
#include "vector3.hpp"
#include "tetrahedron.hpp"
#include "hexahedron.hpp"

void deb(Vector3 vec)
{

}

int main()
{
   Vector3 a(1, 2, 3);
   Vector3 b(3, 4, 5);
   Vector3 c(6, 4, 8);
   Vector3 d(0, 3, 7);
   Vector3 e(4, 6, 3);
   Vector3 f(0, 4, 8);
   Vector3 g(2, 3, 4);
   Vector3 h(1, 5, 3);

   Tetrahedron tetr;

   std::vector<Vector3> vec = { a, b, c, d, e, f, g, h };

   Hexahedron hex(vec);

   vec[0] = Vector3(0, 0, 0);

   auto deb = a + b - c;

}
