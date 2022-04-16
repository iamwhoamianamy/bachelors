#pragma once
#include <ostream>
#include "real.hpp"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

class Complex
{
public:
   real r;
   real i;

   __all__ Complex(real r = 0, real i = 0);
   __all__ Complex(const Complex& complex);

   __all__ Complex& operator+=(const Complex& rhs);
   __all__ Complex& operator*=(const Complex& rhs);

   __all__ Complex& operator=(const Complex& complex);
   __all__ Complex& operator+(const Complex& rhs) const;
   __all__ Complex& operator-(const Complex& rhs) const;
   __all__ Complex& operator*(const Complex& rhs) const;
   __all__ Complex& operator*(const real rhs) const;
};

Complex& operator*(real lhs, const Complex& rhs);
std::ostream& operator<<(std::ostream& os, const Complex& complex);
