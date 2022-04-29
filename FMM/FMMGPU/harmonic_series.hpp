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
   void setHarmonic(int l, int m, T value);
   size_t order() const;
   size_t elemCount() const;

   const T& getReal(int l, int m) const;
   const T& getImag(int l, int m) const;

   static T getRealFactor(int l, int m);
   static T getImagFactor(int l, int m);

   std::vector<T>& data();
   const std::vector<T>& data() const;

   static size_t lmToIndex(int l, int m);

   HarmonicSeries<T>(HarmonicSeries<T>&& harmonicSeries) noexcept;
   HarmonicSeries<T>(const HarmonicSeries<T>& harmonicSeries);
   HarmonicSeries<T>& operator=(HarmonicSeries<T>&& harmonicSeries) noexcept;
   HarmonicSeries<T>& operator=(const HarmonicSeries<T>& harmonicSeries);

   HarmonicSeries<T>(std::vector<T>&& harmonicData) noexcept;
   HarmonicSeries<T>(const std::vector<T>& harmonicData);

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
inline void HarmonicSeries<T>::setHarmonic(int l, int m, T value)
{
   _data[l * l + l + m] = value;
}

template<class T>
inline size_t HarmonicSeries<T>::order() const
{
   return _order;
}

template<class T>
inline size_t HarmonicSeries<T>::elemCount() const
{
   return (_order + 1) * (_order + 1);
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
inline T HarmonicSeries<T>::getRealFactor(int l, int m)
{
   return math::R_SQRT_2 * (m != 0) + (m == 0);
}

template<class T>
inline T HarmonicSeries<T>::getImagFactor(int l, int m)
{
   return math::sign(m) * math::R_SQRT_2;
}

template<class T>
inline std::vector<T>& HarmonicSeries<T>::data()
{
   return _data;
}

template<class T>
inline const std::vector<T>& HarmonicSeries<T>::data() const
{
   return _data;
}

template<class T>
inline size_t HarmonicSeries<T>::lmToIndex(int l, int m)
{
   return l * l + l + m;
}

template<class T>
inline HarmonicSeries<T>& HarmonicSeries<T>::operator=(
   HarmonicSeries<T>&& harmonicSeries) noexcept
{
   _data = std::move(harmonicSeries._data);
   _order = harmonicSeries._order;
   return *this;
}

template<class T>
inline HarmonicSeries<T>& HarmonicSeries<T>::operator=(
   const HarmonicSeries<T>& harmonicSeries)
{
   _data = harmonicSeries._data;
   _order = harmonicSeries._order;
   return *this;
}

template<class T>
inline void HarmonicSeries<T>::add(
   const HarmonicSeries<T>& harmonicSeries)
{
   for(size_t i = 0; i < _data.size(); i++)
   {
      _data[i] += harmonicSeries._data[i];
   }
}

template<class T>
inline void HarmonicSeries<T>::subtract(
   const HarmonicSeries<T>& harmonicSeries)
{
   for(size_t i = 0; i < _data.size(); i++)
   {
      _data[i] -= harmonicSeries._data[i];
   }
}

template<class T>
inline HarmonicSeries<T>::HarmonicSeries(
   HarmonicSeries<T>&& harmonicSeries) noexcept
{
   _data = std::move(harmonicSeries._data);
   _order = harmonicSeries._order;
}

template<class T>
inline HarmonicSeries<T>::HarmonicSeries(
   const HarmonicSeries<T>& harmonicSeries)
{
   _data = harmonicSeries._data;
   _order = harmonicSeries._order;
}

template<class T>
inline HarmonicSeries<T>::HarmonicSeries(
   std::vector<T>&& harmonicData) noexcept
{
   _order = sqrt(harmonicData.size() - 1);
   _data = std::move(harmonicData);
}

template<class T>
inline HarmonicSeries<T>::HarmonicSeries(
   const std::vector<T>& harmonicData)
{
   _order = sqrt(harmonicData.size() - 1);
   _data = harmonicData;
}
