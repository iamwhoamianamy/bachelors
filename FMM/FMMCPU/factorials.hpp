#pragma once
#include <vector>
#include "real.hpp"

class Factorials
{
private:
   std::vector<real> _factorials;
public:
   const size_t maxFactorialNum = 100;
   Factorials();
   size_t operator[](size_t i) const;
};

