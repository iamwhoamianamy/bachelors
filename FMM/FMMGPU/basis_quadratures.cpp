#include <fstream>
#include "basis_quadratures.hpp"
#include "exeptions.hpp"

using namespace std;

BasisQuadratures::BasisQuadratures()
{
   _order = 0;
}

void BasisQuadratures::InitFromTXT(string coordsFileName, string weightsFileName)
{
   InitCoordinatesFromTXT(coordsFileName);
   InitWeightsFromTXT(weightsFileName);
}

int BasisQuadratures::order() const
{
   return _order;
}

const Vector3& BasisQuadratures::values(size_t i) const
{
   return _values[i];
}

const real BasisQuadratures::w(size_t i) const
{
   return _w[i];
}

void BasisQuadratures::InitCoordinatesFromTXT(std::string coordsFileName)
{
   ifstream fin(coordsFileName, ios_base::in);

   if(fin.fail())
      throw Exeption();

   fin.exceptions(ifstream::badbit | ifstream::failbit);

   try
   {
      real x_coord, y_coord, z_coord;

      while(fin >> x_coord && fin >> y_coord && fin >> z_coord)
      {
         _values.push_back(Vector3(x_coord, y_coord, z_coord));
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
         _order = _values.size();
         fin.close();
      }
   }
}

void BasisQuadratures::InitWeightsFromTXT(std::string weightsFileName)
{
   ifstream fin(weightsFileName, ios_base::in);

   if(fin.fail())
      throw Exeption("No such file!");

   fin.exceptions(ifstream::badbit | ifstream::failbit);

   try
   {
      real w_value;

      while(fin >> w_value)
      {
         _w.push_back(w_value);
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
         if(_w.size() != _order)
            throw Exeption("Bad data!");

         fin.close();
      }
   }
}

