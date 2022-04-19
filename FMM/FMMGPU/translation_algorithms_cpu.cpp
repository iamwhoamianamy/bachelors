#include "translation_algorithms.hpp"
#include "kernels.cuh"
#include "cblas.h"

void kernels::translateAllCPU(
   Vector3* result,
   const real* a,
   const Vector3* b,
   size_t harmonicCount, size_t harmonicOrder)
{
   size_t harmonicLength = (harmonicOrder + 1) * (harmonicOrder + 1);

   for(size_t harmonicId = 0; harmonicId < harmonicCount; harmonicId++)
   {
      size_t harmonicBegin = harmonicLength * harmonicId;
      size_t harmonicEnd = harmonicBegin + harmonicLength;

      for(int l = 0; l <= harmonicOrder; l++)
      {
         Vector3 zeroRes = 0;

         for(int lambda = 0; lambda <= l; lambda++)
         {
            int dl = l - lambda;

            for(int mu = -lambda; mu <= lambda; mu++)
               //for(int mu = std::max(-(l - lambda), -lambda); mu <= std::min(l - lambda, lambda); mu++)
            {
               if(-dl <= mu && mu <= +dl)
               {
                  zeroRes += getHarmonic(b, harmonicBegin, lambda, mu) *
                     getHarmonic(a, harmonicBegin, dl, mu) *
                     Harmonics::strangeFactor(0, mu);
               }
            }
         }

         result[lmToIndex(harmonicBegin, l, 0)] = zeroRes;

         for(int m = 1; m <= l; m++)
         {
            Vector3 realRes = 0;
            Vector3 imagRes = 0;

            for(int lambda = 0; lambda <= l; lambda++)
            {
               int dl = l - lambda;

               for(int mu = -lambda; mu <= lambda; mu++)
               {
                  int dm = m - mu;
                  int dnm = -m - mu;

                  Vector3 RR = getReal(b, harmonicBegin, lambda, mu);
                  Vector3 IR = getImag(b, harmonicBegin, lambda, mu);

                  real RM = 0;
                  real IM = 0;

                  real RnM = 0;
                  real InM = 0;

                  if(-dl <= dm && dm <= dl)
                  {
                     RM = getReal(a, harmonicBegin, dl, dm);
                     IM = getImag(a, harmonicBegin, dl, dm);

                     realRes += (RR * RM - IR * IM) * Harmonics::strangeFactor(m, mu);
                     imagRes += (RR * IM + IR * RM) * Harmonics::strangeFactor(m, mu);
                  }

                  if(-dl <= dnm && dnm <= dl)
                  {
                     RnM = getReal(a, harmonicBegin, dl, dnm);
                     InM = getImag(a, harmonicBegin, dl, dnm);

                     realRes += (RR * RnM - IR * InM) * Harmonics::strangeFactor(-m, mu);
                     imagRes -= (RR * InM + IR * RnM) * Harmonics::strangeFactor(-m, mu);
                  }
               }
            }

            result[lmToIndex(harmonicBegin, l, m)] = realRes * math::R_SQRT_2;
            result[lmToIndex(harmonicBegin, l, -m)] = imagRes * math::R_SQRT_2;
         }
      }
   }
}

void kernels::translateAllCPUMatrix(
   std::vector<Complex>& result,
   const std::vector<Complex>& a,
   const std::vector<Complex>& b,
   size_t harmonicCount,
   size_t harmonicOrder)
{
   result = math::multMatricesAsVectors(
      a,
      b,
      (harmonicOrder + 1) * (harmonicOrder + 1),
      harmonicCount,
      (harmonicOrder + 1) * (harmonicOrder + 1));
}


void kernels::translateAllCPUMatrixBLAS(
   Complex* result,
   const Complex* a,
   const Complex* b,
   size_t harmonicCount,
   size_t harmonicOrder)
{
   size_t harmonicLength = (harmonicOrder + 1) * (harmonicOrder + 1);

   size_t harLenPadded = math::nextDevisible(
      harmonicLength,
      THREADS_PER_BLOCK);

   size_t harCountPadded = math::nextDevisible(
      harmonicCount,
      THREADS_PER_BLOCK);

   int m = harLenPadded;
   int k = harLenPadded;
   int n = harCountPadded;
   int lda = m, ldb = k, ldc = m;
   const Complex alf = make_cuComplex(1, 0);
   const Complex bet = make_cuComplex(0, 0);
   const Complex* alpha = &alf;
   const Complex* beta = &bet;

   cblas_cgemm(CBLAS_ORDER::CblasColMajor,
               CBLAS_TRANSPOSE::CblasNoTrans,
               CBLAS_TRANSPOSE::CblasNoTrans,
               m, n, k,
               (float*)alpha, 
               (float*)b, ldb, (float*)a, lda,
               (float*)beta,
               (float*)result, ldc);
}