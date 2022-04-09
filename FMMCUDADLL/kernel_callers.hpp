#pragma once
#include <vector>
#include "../FMM/FMMCPU/real.hpp"

namespace kernels
{
   std::vector<real> __declspec(dllexport)
      addVectors(const std::vector<real>& a,
                 const std::vector<real>& b);
}