#pragma once
#include <vector>
#include "real.hpp"
#include "math.hpp"

template <class T>
struct HarmonicSeries
{
private:
   std::vector<T> _data;
   size_t _n;
public:
   HarmonicSeries();
   HarmonicSeries(int n);

   T& getHarmonic(int l, int m);
   const T& getHarmonic(int l, int m) const;
   size_t size() const;

   const T& getReal(int l, int m) const;
   const T& getImag(int l, int m) const;

   HarmonicSeries<T>(HarmonicSeries<T>&& harmonicSeries) noexcept;
   HarmonicSeries<T>(const HarmonicSeries<T>& harmonicSeries) noexcept;
   HarmonicSeries<T>& operator=(HarmonicSeries<T>&& harmonicSeries) noexcept;
   HarmonicSeries<T>& operator=(const HarmonicSeries<T>& harmonicSeries);

   void add(const HarmonicSeries<T>& harmonicSeries);
   void subtract(const HarmonicSeries<T>& harmonicSeries);
};

template<class T>
inline HarmonicSeries<T>::HarmonicSeries()
{
}

template<class T>
inline HarmonicSeries<T>::HarmonicSeries(int n)
{
   _n = n;
   _data = std::vector<T>((n + 1) * (n + 1));
}

template<class T>
inline T& HarmonicSeries<T>::getHarmonic(int l, int m)
{
   return _data[l * l + l + m];
}

template<class T>
inline const T& HarmonicSeries<T>::getHarmonic(int l, int m) const
{
   return _data[l * l + l + m];
}

template<class T>
inline size_t HarmonicSeries<T>::size() const
{
   return _n;
}

template<class T>
inline const T& HarmonicSeries<T>::getReal(int l, int m) const
{
   return getHarmonic(l, abs(m)) * 
      (math::R_SQRT_2 * (m != 0) + (m == 0));
}

template<class T>
inline const T& HarmonicSeries<T>::getImag(int l, int m) const
{
   return math::sign(m) * getHarmonic(l, -abs(m)) * math::R_SQRT_2 * (m != 0);
}

template<class T>
inline HarmonicSeries<T>& HarmonicSeries<T>::operator=(HarmonicSeries<T>&& harmonicSeries) noexcept
{
   _data = std::move(harmonicSeries._data);
   _n = harmonicSeries._n;
   return *this;
}

template<class T>
inline HarmonicSeries<T>& HarmonicSeries<T>::operator=(const HarmonicSeries<T>& harmonicSeries)
{
   _data = harmonicSeries._data;
   _n = harmonicSeries._n;
   return *this;
}

template<class T>
inline void HarmonicSeries<T>::add(const HarmonicSeries<T>& harmonicSeries)
{
   for(size_t i = 0; i < _data.size(); i++)
   {
      _data[i] += harmonicSeries._data[i];
   }
}

template<class T>
inline void HarmonicSeries<T>::subtract(const HarmonicSeries<T>& harmonicSeries)
{
   for(size_t i = 0; i < _data.size(); i++)
   {
      _data[i] -= harmonicSeries._data[i];
   }
}

template<class T>
inline HarmonicSeries<T>::HarmonicSeries(HarmonicSeries<T>&& harmonicSeries) noexcept
{
   _data = std::move(harmonicSeries._data);
   _n = harmonicSeries._n;
}

template<class T>
inline HarmonicSeries<T>::HarmonicSeries(const HarmonicSeries<T>& harmonicSeries) noexcept
{
   _data = harmonicSeries._data;
   _n = harmonicSeries._n;
}
