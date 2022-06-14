#include "blass_callers.hpp"

void blas::multMatrices(
   const CBLAS_ORDER order,
   const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB,
   const blasint m, const blasint n, const blasint k,
   const real alpha,
   const real* a, const blasint lda,
   const real* b, const blasint ldb,
   const real beta,
   real* c, const blasint ldc)
{
#ifdef REAL_IS_FLOAT

   cblas_sgemm(
      order, transA, transB, m, n, k, alpha, 
      a, lda, b, ldb, beta, c, ldc);

#endif

#ifdef REAL_IS_DOUBLE
   
   cblas_dgemm(
      order, transA, transB, m, n, k, alpha,
      a, lda, b, ldb, beta, c, ldc);

#endif
}

void blas::multComplexMatrices(
   CBLAS_ORDER order,
   CBLAS_TRANSPOSE transA,
   CBLAS_TRANSPOSE transB,
   blasint m, blasint n, blasint k,
   real* alpha,
   real* a, blasint lda,
   real* b, blasint ldb,
   real* beta,
   real* c, blasint ldc)
{
#ifdef REAL_IS_FLOAT

   cblas_cgemm(
      order, transA, transB, m, n, k, alpha,
      a, lda, b, ldb, beta, c, ldc);

#endif

#ifdef REAL_IS_DOUBLE

   cblas_zgemm(
      order, transA, transB, m, n, k, alpha,
      a, lda, b, ldb, beta, c, ldc);

#endif
}

void blas::multMatricesCUBLAS(
   cublasHandle_t handle,
   cublasOperation_t transa,
   cublasOperation_t transb,
   int m, int n, int k,
   const real* alpha,
   const real* a, int lda,
   const real* b, int ldb,
   const real* beta,
   real* c, int ldc)
{
#ifdef REAL_IS_FLOAT

   cublasSgemm_v2(
      handle, transa, transb, m, n, k, alpha,
      a, lda, b, ldb, beta, c, ldc);

#endif

#ifdef REAL_IS_DOUBLE

   cublasDgemm_v2(
      handle, transa, transb, m, n, k, alpha,
      a, lda, b, ldb, beta, c, ldc);

#endif
}

void blas::multMatricesStridedBatchedCUBLAS(
   cublasHandle_t handle,
   cublasOperation_t transa,
   cublasOperation_t transb,
   int m, int n, int k,
   const real* alpha,
   const real* A,int lda, long long strideA,
   const real* B, int ldb, long long strideB,
   const real* beta,
   real* C, int ldc, long long strideC,
   int batchCount)
{
#ifdef REAL_IS_FLOAT

   cublasSgemmStridedBatched(
      handle, transa, transb, m, n, k, alpha,
      A, lda, strideA, B, ldb, strideB, beta,
      C, ldc, strideC, batchCount);

#endif

#ifdef REAL_IS_DOUBLE

   cublasDgemmStridedBatched(
      handle, transa, transb, m, n, k, alpha,
      A, lda, strideA, B, ldb, strideB, beta,
      C, ldc, strideC, batchCount);

#endif
}

void blas::copyVector(
   blasint n,
   const real* x,
   blasint incx,
   real* y,
   blasint incy)
{
#ifdef REAL_IS_FLOAT

   cblas_scopy(n, x, incx, y, incy);

#endif

#ifdef REAL_IS_DOUBLE

   cblas_dcopy(n, x, incx, y, incy);

#endif
}

void blas::addVectorToVector(
   blasint n,
   real alpha,
   const real* x,
   blasint incx,
   real* y,
   blasint incy)
{
#ifdef REAL_IS_FLOAT

   cblas_saxpy(n, alpha, x, incx, y, incy);

#endif

#ifdef REAL_IS_DOUBLE

   cblas_daxpy(n, alpha, x, incx, y, incy);

#endif
}
