#pragma once
#include <vector>

template <class T>
struct HarmonicSeries
{
   std::vector<std::vector<T>> data;
   HarmonicSeries();
   HarmonicSeries(int n);
   T& getHarmonic(int l, int m);
   std::vector<T>& operator[](int i);
   HarmonicSeries<T>(HarmonicSeries<T>&& harmonicSeries) noexcept;
   HarmonicSeries<T>(const HarmonicSeries<T>& harmonicSeries) noexcept;
   HarmonicSeries<T>& operator=(HarmonicSeries<T>&& harmonicSeries) noexcept;
   HarmonicSeries<T>& operator=(const HarmonicSeries<T>& harmonicSeries);
};

template<class T>
inline HarmonicSeries<T>::HarmonicSeries()
{
}

template<class T>
inline HarmonicSeries<T>::HarmonicSeries(int n)
{
   data = std::vector<std::vector<T>>(n);

   for(size_t i = 0; i < n; i++)
   {
      data[i] = std::vector<T>(2 * i + 1);
   }
}

template<class T>
inline T& HarmonicSeries<T>::getHarmonic(int l, int m)
{
   return data[l][data[l].size() / 2 + m];
}

template<class T>
inline std::vector<T>& HarmonicSeries<T>::operator[](int i)
{
   return data[i];
}

template<class T>
inline HarmonicSeries<T>& HarmonicSeries<T>::operator=(HarmonicSeries<T>&& harmonicSeries) noexcept
{
   data = std::move(harmonicSeries.data);
   return *this;
}

template<class T>
inline HarmonicSeries<T>& HarmonicSeries<T>::operator=(const HarmonicSeries<T>& harmonicSeries)
{
   data = harmonicSeries.data;
   return *this;
}

template<class T>
inline HarmonicSeries<T>::HarmonicSeries(HarmonicSeries<T>&& harmonicSeries) noexcept
{
   data = std::move(harmonicSeries.data);
}

template<class T>
inline HarmonicSeries<T>::HarmonicSeries(const HarmonicSeries<T>& harmonicSeries) noexcept
{
   data = harmonicSeries.data;
}
