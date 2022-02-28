#include <fstream>
#include "basis_quadratures.hpp"
#include "exeptions.hpp"

using namespace std;

BasisQuadratures::BasisQuadratures()
{
   order = 0;
}

void BasisQuadratures::InitFromTXT(string coordsFileName, string weightsFileName)
{
   InitCoordinatesFromTXT(coordsFileName);
   InitWeightsFromTXT(weightsFileName);
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
         x.push_back(x_coord);
         y.push_back(y_coord);
         z.push_back(z_coord);
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
         order = x.size();
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
         w.push_back(w_value);
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
         if(w.size() != order)
            throw Exeption("Bad data!");

         fin.close();
      }
   }
}
