#include "exeptions.h"

namespace trq = triangle_quadratures;

std::ostream& trq::operator<<(std::ostream& out, const trq::Exeption& exeption)
{
   out << exeption.message;
   return out;
}

triangle_quadratures::ParsingExeption::ParsingExeption() : Exeption("Parsing exeption!")
{
 
}
