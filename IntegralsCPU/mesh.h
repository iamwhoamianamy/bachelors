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
      Mesh();
      //void InitFromTXT(string fileName);
      void InitFromOBJ(string fileName);
      //Mesh(Mesh&& mesh) noexcept;

      Triangle GetTriangle(int index) const;
      int TriangleCount() const;

      //Mesh& operator=(Mesh&& mesh) noexcept;
   };
}
