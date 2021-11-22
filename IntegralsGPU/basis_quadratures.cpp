#include <fstream>
#include "basis_quadratures.h"
#include "exeptions.h"

using namespace std;

namespace triangle_quadratures
{
   BasisQuadratures::BasisQuadratures()
   {
      order = 0;
   }

   BasisQuadratures::BasisQuadratures(int order) : order(order)
   {
      switch(order)
      {
         case 0:
         {
            break;
         }
         case 6:
         {
            x = { 0.659027622374092, 0.659027622374092, 0.231933368553031, 0.231933368553031, 0.109039009072877, 0.109039009072877 };
            y = { 0.231933368553031, 0.109039009072877, 0.659027622374092, 0.109039009072877, 0.659027622374092, 0.231933368553031 };
            w = { 0.166666666666667, 0.166666666666667, 0.166666666666667, 0.166666666666667, 0.166666666666667, 0.166666666666667 };
            
            break;
         }
         default:
         {
            throw NotImplementedExeption();
         }
      }
   }

   void BasisQuadratures::InitFromTXT(string coordsFileName, string weightsFileName)
   {
      ifstream fin;
      fin.open(coordsFileName, ios_base::in);

      if(fin.fail())
         throw Exeption("No such file!");

      fin.exceptions(ifstream::badbit | ifstream::failbit);

      try
      {
         real x_coord, y_coord;

         while(fin >> x_coord && fin >> y_coord)
         {
            x.push_back(x_coord);
            y.push_back(y_coord);
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

      fin.open(weightsFileName, ios_base::in);

      if(fin.fail())
         throw Exeption("No such file!");

      fin.exceptions(ifstream::badbit | ifstream::failbit);

      try
      {
         real w_value, y_coord;

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
};