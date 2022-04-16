#include "complex.cuh"
namespace math
{
   Complex::Complex() :
      _r(0), _i(0)
   {
   }

   Complex::Complex(real value) :
      _r(value), _i(0)
   {
   }

   Complex::Complex(real r, real i) :
      _r(r), _i(i)
   {

   }

   Complex::Complex(const Complex& complex) :
      _r(complex.r()), _i(complex.i())
   {
   }

   Complex::Complex(Complex&& complex) noexcept
   {
      _r = std::move(complex._r);
      _i = std::move(complex._i);
   }

   real Complex::r() const
   {
      return _r;
   }

   real Complex::i() const
   {
      return _i;
   }

   Complex& Complex::operator+=(const Complex& rhs)
   {
      _r += rhs.r();
      _i += rhs.i();

      return *this;
   }

   Complex& Complex::operator*=(const Complex& rhs)
   {
      real newR = _r * rhs.r() - _i * rhs.i();
      real newI = _i * rhs.r() + _r * rhs.i();

      _r = newR;
      _i = newI;

      return *this;
   }

   Complex& Complex::operator=(const Complex& complex)
   {
      _r = complex.r();
      _i = complex.i();

      return *this;
   }

   Complex& Complex::operator=(Complex&& complex) noexcept
   {
      _r = std::move(complex._r);
      _i = std::move(complex._i);

      return *this;
   }

   Complex& Complex::operator+(const Complex& rhs) const
   {
      return Complex(_r + rhs.r(), _i + rhs.i());
   }

   Complex& Complex::operator-(const Complex& rhs) const
   {
      return Complex(_r - rhs.r(), _i - rhs.i());
   }

   Complex& Complex::operator*(const Complex& rhs) const
   {
      return Complex(_r * rhs.r() - _i * rhs.i(), _i * rhs.r() + _r * rhs.i());
   }

   Complex& Complex::operator*(const real rhs) const
   {
      return Complex(_r * rhs, _i * rhs);
   }
}

std::ostream& operator<<(std::ostream& os, const math::Complex& complex)
{
   os << "(" << complex.r() << "," << complex.i() << ")";
   return os;
}