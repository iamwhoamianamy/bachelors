#pragma once
#include "complex.cuh"

Complex::Complex(real r, real i) :
   r(r), i(i)
{

}

Complex::Complex(const Complex& complex) :
   r(complex.r), i(complex.i)
{
}

__all__ Complex& Complex::operator+=(const Complex& rhs)
{
   r += rhs.r;
   i += rhs.i;

   return *this;
}

__all__ Complex& Complex::operator*=(const Complex& rhs)
{
   r = r * rhs.r - i * rhs.i;
   i = i * rhs.r + r * rhs.i;

   return *this;
}

Complex& Complex::operator=(const Complex& complex)
{
   r = complex.r;
   i = complex.i;

   return *this;
}

Complex& Complex::operator+(const Complex& rhs) const
{
   return Complex(r + rhs.r, i + rhs.i);
}

Complex& Complex::operator-(const Complex& rhs) const
{
   return Complex(r - rhs.r, i - rhs.i);
}

Complex& Complex::operator*(const Complex& rhs) const
{
   return Complex(r * rhs.r - i * rhs.i, i * rhs.r + r * rhs.i);
}

Complex& Complex::operator*(const real rhs) const
{
   return Complex(r * rhs, i * rhs);
}

Complex& operator*(real lhs, const Complex& rhs)
{
   return rhs * lhs;
}

std::ostream& operator<<(std::ostream& os, const Complex& complex)
{
   os << "(" << complex.r << "," << complex.i << ")";
   return os;
}
