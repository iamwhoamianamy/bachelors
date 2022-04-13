#pragma once
#include <vector>
#include "real.hpp"
#include "math.hpp"

template <class T>
struct HarmonicSeries
{
private:
   std::vector<T> _data;
   size_t _order;
public:
   HarmonicSeries();
   HarmonicSeries(int order);

   T& getHarmonic(int order);
   const T& getHarmonic(int order) const;

   T& getHarmonic(int l, int m);
   const T& getHarmonic(int l, int m) const;
   size_t order() const;

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
inline HarmonicSeries<T>::HarmonicSeries(int order)
{
   _order = order;
   _data = std::vector<T>((order + 1) * (order + 1));
}

template<class T>
inline T& HarmonicSeries<T>::getHarmonic(int order)
{
   return _data[order];
}

template<class T>
inline const T& HarmonicSeries<T>::getHarmonic(int order) const
{
   return _data[order];
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
inline size_t HarmonicSeries<T>::order() const
{
   return _order;
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
   _order = harmonicSeries._order;
   return *this;
}

template<class T>
inline HarmonicSeries<T>& HarmonicSeries<T>::operator=(const HarmonicSeries<T>& harmonicSeries)
{
   _data = harmonicSeries._data;
   _order = harmonicSeries._order;
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
   _order = harmonicSeries._order;
}

template<class T>
inline HarmonicSeries<T>::HarmonicSeries(const HarmonicSeries<T>& harmonicSeries) noexcept
{
   _data = harmonicSeries._data;
   _order = harmonicSeries._order;
}
