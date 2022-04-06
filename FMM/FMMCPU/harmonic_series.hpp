#pragma once
#include <vector>
#include "real.hpp"
#include "math.hpp"

template <class T>
struct HarmonicSeries
{
private:
   std::vector<std::vector<T>> _data;
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
   _data = std::vector<std::vector<T>>(n);

   for(size_t i = 0; i < n; i++)
   {
      _data[i] = std::vector<T>(2 * i + 1);
   }
}

template<class T>
inline T& HarmonicSeries<T>::getHarmonic(int l, int m)
{
   return _data[l][l + m];
}

template<class T>
inline const T& HarmonicSeries<T>::getHarmonic(int l, int m) const
{
   return _data[l][l + m];
}

template<class T>
inline size_t HarmonicSeries<T>::size() const
{
   return _data.size();
}

template<class T>
inline const T& HarmonicSeries<T>::getReal(int l, int m) const
{
   if(m == 0) return getHarmonic(l, 0);
   return getHarmonic(l, abs(m)) * math::R_SQRT_2;
}

template<class T>
inline const T& HarmonicSeries<T>::getImag(int l, int m) const
{
   if(m == 0) return T(0);
   return math::sign(m) * getHarmonic(l, -abs(m)) * math::R_SQRT_2;
}

template<class T>
inline HarmonicSeries<T>& HarmonicSeries<T>::operator=(HarmonicSeries<T>&& harmonicSeries) noexcept
{
   _data = std::move(harmonicSeries._data);
   return *this;
}

template<class T>
inline HarmonicSeries<T>& HarmonicSeries<T>::operator=(const HarmonicSeries<T>& harmonicSeries)
{
   _data = harmonicSeries._data;
   return *this;
}

template<class T>
inline void HarmonicSeries<T>::add(const HarmonicSeries<T>& harmonicSeries)
{
   for(int l = 0; l < _data.size(); l++)
   {
      for(int m = -l; m <= l; m++)
      {
         getHarmonic(l, m) += harmonicSeries.getHarmonic(l, m);
      }
   }
}

template<class T>
inline void HarmonicSeries<T>::subtract(const HarmonicSeries<T>& harmonicSeries)
{
   for(int l = 0; l < _data.size(); l++)
   {
      for(int m = -l; m <= l; m++)
      {
         getHarmonic(l, m) -= harmonicSeries.getHarmonic(l, m);
      }
   }
}

template<class T>
inline HarmonicSeries<T>::HarmonicSeries(HarmonicSeries<T>&& harmonicSeries) noexcept
{
   _data = std::move(harmonicSeries._data);
}

template<class T>
inline HarmonicSeries<T>::HarmonicSeries(const HarmonicSeries<T>& harmonicSeries) noexcept
{
   _data = harmonicSeries._data;
}
