#pragma once
#include <vector>

class Factorials
{
private:
   std::vector<size_t> _factorials;
public:
   const size_t maxFactorialNum = 100;
   Factorials();
   size_t operator[](size_t i) const;
};

