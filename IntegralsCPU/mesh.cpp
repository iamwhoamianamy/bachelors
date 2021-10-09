#include <fstream>
#include "mesh.h"
#include "exeptions.h"

namespace triangle_quadratures
{
   Mesh::Mesh() : _triangles_count(0)
   {

   }

   /*void Mesh::InitFromTXT(string fileName)
   {
      ifstream fin;
      fin.open(fileName, ios_base::in);

      if(fin.fail())
         throw FileExeption();

      fin >> _triangles_count;
      _vertices = vector<Vector3>(_triangles_count);

      for(size_t i = 0; i < _triangles_count; i++)
      {
         double x, y, z;
         fin >> x >> y >> z;
         _vertices[i] = Vector3(x, y, z);
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
   }*/

   void Mesh::InitFromOBJ(string fileName)
   {
      ifstream fin;
      fin.open(fileName, ios_base::in);

      if(fin.fail())
         throw Exeption("No such file!");

      fin.exceptions(ifstream::badbit | ifstream::failbit);

      try
      {
         string letter;

         while(fin >> letter)
         {

            if(letter == "v")
            {
               double x, y, z;
               fin >> x >> y >> z;

               _vertices.push_back(Vector3(x, y, z));
            }
            if(letter == "f")
            {
               vector<int> temp(3);

               int i0, i1, i2;
               fin >> i0 >> i1 >> i2;

               temp[0] = i0 - 1;
               temp[1] = i1 - 1;
               temp[2] = i2 - 1;

               if(temp[0] < 0 || temp[1] < 0 || temp[2] < 0 ||
                  temp[0] >= _vertices.size() || temp[1] >= _vertices.size() || temp[2] >= _vertices.size())
                  throw Exeption("OBJ file is damaged!");

               _faces_ind.push_back(temp);
            }
         }

      }
      catch(ifstream::failure e)
      {
         if(!fin.eof())
         {
            fin.close();
            throw ParsingExeption();
         }
         else
         {
            _triangles_count = _faces_ind.size();
            fin.close();
         }
      }
   }

   //Mesh::Mesh(Mesh&& mesh) noexcept
   //{
   //   _triangles_count = mesh._triangles_count;
   //   _vertices = std::move(mesh._vertices);
   //   _faces_ind = std::move(mesh._faces_ind);
   //}

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

   //Mesh& Mesh::operator=(Mesh&& mesh) noexcept
   //{
   //   _triangles_count = mesh._triangles_count;
   //   _vertices = std::move(mesh._vertices);
   //   _faces_ind = std::move(mesh._faces_ind);
   //   return *this;
   //}
}
