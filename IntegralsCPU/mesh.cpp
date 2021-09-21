#include <fstream>
#include "mesh.h"
#include "exeptions.h"

namespace triangle_quadratures
{
   Mesh::Mesh(string fileName)
   {
      ifstream fin;
      fin.open(fileName, ios_base::in);

      if(fin.fail())
         throw FileExeption();

      fin >> _triangles_count;
      _vertices = vector<Point>(_triangles_count);

      for(size_t i = 0; i < _triangles_count; i++)
      {
         double x, y, z;
         fin >> x >> y >> z;
         _vertices[i] = Point(x, y, z);
      }

      int faces_count;
      fin >> faces_count;
      _faces_ind = vector<vector<int>>(faces_count, vector<int>(3));

      for(size_t i = 0; i < faces_count; i++)
      {
         int i0, i1, i2;
         fin >> i0 >> i1 >> i2;

         _faces_ind[i][0] = i0 - 1;
         _faces_ind[i][1] = i1 - 1;
         _faces_ind[i][2] = i2 - 1;
      }

      fin.close();
   }

   Triangle Mesh::GetTriangle(int index) const
   {
      if(index >= _triangles_count)
         throw RangeExeption();

      return Triangle(_vertices[_faces_ind[index][0]],
                      _vertices[_faces_ind[index][1]],
                      _vertices[_faces_ind[index][2]]);
   }
   int Mesh::TriangleCount() const
   {
      return _triangles_count;
   }
}
