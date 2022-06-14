#pragma once
#include <vector>
#include <ostream>
#include "real.hpp"
#include "cuda_runtime.h"
#include "typedefs.hpp"
#include "vector3.cuh"
#include "box.hpp"

namespace math
{
   constexpr real PI = 3.14159265359;
   constexpr real MU0 = 1.2566370614e-6;
   constexpr real eps = 1e-3;
   const real SQRT_2 = sqrt(2.0);
   __device__ constexpr real R_SQRT_2 = 0.70710678118;

   Vector3 cylindricToCartesian(const Vector3& point);

   real calcFactorial(int n);
   real calcBinomial(int k, int n);
   real randBetween(real min, real max);
   size_t nextDevisible(const size_t number, const size_t devidor);
   Box getBoundingBox(const std::vector<Vector3>& points);

   template <class T>
   constexpr __all__ int sign(T val)
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
   std::vector<T> mult(
      const std::vector<T>& a,
      const Matrix<T>& b)
   {
      std::vector<T> res(b[0].size());

      for(size_t column = 0; column < b[0].size(); column++)
      {
         for(size_t row = 0; row < a.size(); row++)
         {
            res[column] += a[row] * b[row][column];
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
               res[i * bWidth + j] += a[i * aWidth + k] * b[k * bWidth + j];
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
      const size_t width,
      const size_t height,
      const size_t padding,
      const int idx)
   {
      const size_t currentWidth = nextDevisible(width, padding);
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
      const size_t width,
      const size_t padding,
      const int idx)
   {
      const size_t currentWidth = nextDevisible(width, padding);
      std::vector<T> res(matrix.begin() + idx * currentWidth,
          matrix.begin() + idx * currentWidth + width);

      return res;
   }

   template<class T>
   std::vector<T> matrixToVector(const Matrix<T>& matrix, const size_t padding)
   {
      const size_t height = math::nextDevisible(matrix.size(), padding);
      const size_t width = math::nextDevisible(matrix[0].size(), padding);
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
      const size_t height,
      const size_t width,
      const size_t padding)
   {
      const size_t currentWidth = math::nextDevisible(width, padding);
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
   
   real max(const std::vector<Vector3>& vec, size_t axis);
   real min(const std::vector<Vector3>& vec, size_t axis);

   template<class T>
   T sum(const std::vector<T>& vec)
   {
      T res = 0;

      for(size_t i = 0; i < vec.size(); i++)
      {
         res += vec[i];
      }

      return res;
   }

   template<class T>
   T mean(const std::vector<T>& vec)
   {
      return sum(vec) / vec.size();
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

__all__ real getReal(const Complex& complex);
__all__ real getImag(const Complex& complex);
__all__ Complex makeComplex(real realPart, real imagPart);
__all__ Complex operator*(const Complex& lhs, const Complex& rhs);
__all__ Complex operator*(const Complex& lhs, real rhs);
__all__ Complex operator+(const Complex& lhs, const Complex& rhs);
__all__ Complex operator-(const Complex& lhs, const Complex& rhs);
__all__ Complex& operator*=(Complex& lhs, const Complex& rhs);
__all__ Complex& operator+=(Complex& lhs, const Complex& rhs);

std::ostream& operator<<(std::ostream& os, const Complex& val);