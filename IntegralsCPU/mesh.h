#pragma once
#include <vector>
#include <string>
#include "triangle.h"

using namespace std;

namespace triangle_quadratures
{
   class Mesh
   {
   private:
      int _triangles_count;
      vector<Point> _vertices;
      vector<vector<int>> _faces_ind;

   public:
      Mesh(string fileName);
      Triangle GetTriangle(int index) const;
      int TriangleCount() const;
   };
}
