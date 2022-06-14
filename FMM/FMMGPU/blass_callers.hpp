#pragma once
#include <cblas.h>
#include "real.hpp"
#include "cublasLt.h"

namespace blas
{
   void multMatrices(
      CBLAS_ORDER order,
      CBLAS_TRANSPOSE transA,
      CBLAS_TRANSPOSE transB,
      blasint m, blasint n, blasint k,
      real alpha,
      const real* a, blasint lda,
      const real* b, blasint ldb,
      real beta,
      real* c, blasint ldc);

   void multMatricesCUBLAS(
      cublasHandle_t handle,
      cublasOperation_t transa,
      cublasOperation_t transb,
      int m, int n, int k,
      const real* alpha,
      const real* a, int lda,
      const real* b, int ldb,
      const real* beta,
      real* c, int ldc);

   void multMatricesStridedBatchedCUBLAS(
      cublasHandle_t handle,
      cublasOperation_t transa,
      cublasOperation_t transb,
      int m, int n, int k,
      const real* alpha,
      const real* A, int lda, long long int strideA,
      const real* B, int ldb, long long int strideB,
      const real* beta, real* C, int ldc, long long int strideC,
      int batchCount);

   void copyVector(
      blasint n,
      const real* x,
      blasint incx,
      real* y,
      blasint incy);

   void addVectorToVector(
      blasint n,
      real alpha,
      const real* x,
      blasint incx, real* y,
      blasint incy);
}
