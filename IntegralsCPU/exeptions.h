#pragma once
#include <string>
#include <ostream>

namespace triangle_quadratures
{
   class NotImplementedExeption{};
   class RangeExeption{};

   class Exeption
   {
   public:
      std::string message;
      Exeption(std::string message) : message(message) {}
   };

   class ParsingExeption : public Exeption
   {
   public:
      ParsingExeption();
   };

   std::ostream& operator <<(std::ostream& out, const Exeption& exeption);
}