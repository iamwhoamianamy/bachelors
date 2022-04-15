#pragma once
#include <vector>
#include "real.hpp"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "typedefs.hpp"

namespace math
{
   const real PI = 3.14159265359;
   const real mu0 = 1.2566370614e-6;
   const real SQRT_2 = sqrt(2.0);
   __device__ const real R_SQRT_2 = 0.70710678118;

   real calcFactorial(int n);
   real calcBinomial(int k, int n);

   template <class T>
   __all__ int sign(T val)
   {
      return (T(0) < val) - (val < T(0));
   }

   template <class T>
   T mult(
      const std::vector<T>& a,
      const std::vector<T>& b)
   {
      T res = 0;

      for(size_t i = 0; i < a.size(); i++)
      {
         res += a[i] * b[i];
      }

      return res;
   }

   template <class T>
   std::vector<T> mult(
      const Matrix<T>& a,
      const std::vector<T>& b)
   {
      std::vector<T> res(b.size());

      for(size_t y = 0; y < a.size(); y++)
      {
         for(size_t x = 0; x < b.size(); x++)
         {
            res[y] += a[y][x] * b[x];
         }
      }

      return res;
   }

   template <class T>
   Matrix<T> mult(
      const Matrix<T>& a,
      const Matrix<T>& b)
   {
      Matrix<T> res(a.size(), std::vector<T>(b[0].size()));

      for(size_t i = 0; i < a.size(); i++)
      {
         for(size_t j = 0; j < b[0].size(); j++)
         {
            for(size_t k = 0; k < b.size(); k++)
            {
               res[i][j] += a[i][k] * b[k][j];
            }
         }
      }

      return res;
   }

   template<class T>
   std::vector<T> getColumn(const Matrix<T>& matrix, int idx)
   {
      std::vector<T> res(matrix.size());

      for(size_t i = 0; i < matrix.size(); i++)
      {
         res[i] = matrix[i][idx];
      }

      return res;
   }

   template <class T>
   std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
   {
      for(size_t i = 0; i < vec.size(); i++)
      {
         os << vec[i] << " ";
      }

      return os;
   }

   template <class T>
   std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix)
   {
      for(size_t i = 0; i < matrix.size(); i++)
      {
         os << matrix[i] << std::endl;
      }

      return os;
   }
}