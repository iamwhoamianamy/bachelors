#pragma once
#include <ostream>
#include "real.hpp"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

class Complex
{
private:
   real _r;
   real _i;

public:
   __all__ Complex();
   __all__ Complex(real value);
   __all__ Complex(real r, real i);
   __all__ Complex(const Complex& complex);
   __all__ Complex(Complex&& complex) noexcept;

   __all__ real r() const;
   __all__ real i() const;

   __all__ Complex& operator+=(const Complex& rhs);
   __all__ Complex& operator*=(const Complex& rhs);

   __all__ Complex& operator=(const Complex& complex);
   __all__ Complex& operator=(Complex&& complex) noexcept;

   __all__ Complex& operator+(const Complex& rhs) const;
   __all__ Complex& operator-(const Complex& rhs) const;
   __all__ Complex& operator*(const Complex& rhs) const;
   __all__ Complex& operator*(const real rhs) const;
};

std::ostream& operator<<(std::ostream& os, const Complex& complex);
