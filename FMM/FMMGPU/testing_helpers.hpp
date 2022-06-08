#pragma once
#include <chrono>

#include "torus.hpp"
#include "basis_quadratures.hpp"
#include "box.hpp"
#include "cylinder.hpp"

namespace test
{
   Torus createTorus();
   Cylinder createCylinder();
   double getTime(void (*f)());
   double getTime(const std::chrono::steady_clock::time_point& start,
                  const std::chrono::steady_clock::time_point& stop);
   BasisQuadratures readTetrahedronBasisQuadratures();
   BasisQuadratures readTriangleBasisQuadratures();
   std::vector<Vector3> createPoints(const Vector3& begin, const Vector3& end, int n);
   std::vector<Vector3> createRandomPoints(const Box& box, int n);
   void printSeparateLine(std::ostream& os, size_t count);
}
