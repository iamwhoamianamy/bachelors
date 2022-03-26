#include "factorials.hpp"

Factorials::Factorials()
{
   _factorials.resize(maxFactorialNum);
   _factorials[0] = 1;

   for(size_t i = 1; i < maxFactorialNum; i++)
   {
      _factorials[i] = _factorials[i - 1] * i;
   }
}

size_t Factorials::operator[](size_t i) const
{
   return _factorials[i];
}
