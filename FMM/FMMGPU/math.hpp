#pragma once
#include <vector>
#include <ostream>
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

   size_t nextDevisible(const size_t number, const size_t devidor);

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
   std::vector<T> multMatricesAsVectors(
      const std::vector<T>& a,
      const std::vector<T>& b,
      size_t aWidth,
      size_t aHeight,
      size_t bWidth)
   {
      std::vector<T> res(aHeight * bWidth);

      for(size_t i = 0; i < aHeight; i++)
      {
         for(size_t j = 0; j < bWidth; j++)
         {
            for(size_t k = 0; k < aWidth; k++)
            {
               res[i * aWidth + j] += a[i * aWidth + k] * b[k * aWidth + j];
            }
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

   template<class T>
   std::vector<T> getColumn(
      const std::vector<T>& matrix,
      size_t width,
      size_t height,
      size_t padding,
      int idx)
   {
      size_t currentWidth = nextDevisible(width, padding);
      std::vector<T> res(height);

      for(size_t i = 0; i < height; i++)
      {
         res[i] = matrix[i * currentWidth + idx];
      }

      return res;
   }

   template<class T>
   std::vector<T> getRow(
      const std::vector<T>& matrix,
      size_t width,
      size_t padding,
      int idx)
   {
      size_t currentWidth = nextDevisible(width, padding);
      std::vector<T> res(matrix.begin() + idx * currentWidth,
          matrix.begin() + idx * currentWidth + width);

      return res;
   }

   template<class T>
   std::vector<T> matrixToVector(const Matrix<T>& matrix, size_t padding)
   {
      size_t height = math::nextDevisible(matrix.size(), padding);
      size_t width = math::nextDevisible(matrix[0].size(), padding);
      std::vector<T> res(height * width);

      for(size_t i = 0; i < matrix.size(); i++)
      {
         for(size_t j = 0; j < matrix[0].size(); j++)
         {
            res[i * width + j] = matrix[i][j];
         }
      }

      return res;
   }

   template<class T>
   Matrix<T> vectorToMatrix(
      const std::vector<T>& vector,
      size_t height,
      size_t width,
      size_t padding)
   {
      size_t currentWidth = math::nextDevisible(width, padding);
      Matrix<T> res(height, std::vector<T>(width));

      for(size_t i = 0; i < height; i++)
      {
         for(size_t j = 0; j < width; j++)
         {
            res[i][j] = vector[i * currentWidth + j];
         }
      }

      return res;
   }

   template<class T>
   void translateSquareMatrix(Matrix<T>& matrix)
   {
      for(size_t i = 0; i < matrix.size(); i++)
      {
         for(size_t j = 0; j < i; j++)
         {
            std::swap(matrix[i][j], matrix[j][i]);
         }
      }
   }
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