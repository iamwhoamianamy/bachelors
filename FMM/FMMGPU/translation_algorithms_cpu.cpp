#include "translation_algorithms.hpp"
#include "kernels.cuh"
#include "cblas.h"
#include "multipole_translator.hpp"
#include "blass_callers.hpp"

#include <omp.h>

void kernels::translateAllCPU(
   Vector3* result,
   const real* a,
   const Vector3* b,
   size_t harmonicCount, size_t harmonicOrder)
{
   size_t harmonicLength = (harmonicOrder + 1) * (harmonicOrder + 1);

#pragma omp parallel for
   for(int harmonicId = 0; harmonicId < harmonicCount; harmonicId++)
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
            {
               if(-dl <= mu && mu <= +dl)
               {
                  zeroRes += getHarmonic(b, harmonicBegin, lambda, mu) *
                     getHarmonic(a, harmonicBegin, dl, mu) *
                     MultipoleTranslator::multipoleTranslationFactor(0, mu);
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

                     realRes += (RR * RM - IR * IM) * 
                        MultipoleTranslator::multipoleTranslationFactor(m, mu);
                     imagRes += (RR * IM + IR * RM) * 
                        MultipoleTranslator::multipoleTranslationFactor(m, mu);
                  }

                  if(-dl <= dnm && dnm <= dl)
                  {
                     RnM = getReal(a, harmonicBegin, dl, dnm);
                     InM = getImag(a, harmonicBegin, dl, dnm);

                     realRes += (RR * RnM - IR * InM) * 
                        MultipoleTranslator::multipoleTranslationFactor(-m, mu);
                     imagRes -= (RR * InM + IR * RnM) * 
                        MultipoleTranslator::multipoleTranslationFactor(-m, mu);
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
   real* result,
   const real* a,
   const real* b,
   size_t harmonicCount,
   size_t harmonicOrder)
{
   size_t harmonicLength = (harmonicOrder + 1) * (harmonicOrder + 1);

   int m = harmonicLength;
   int k = harmonicLength;
   int n = harmonicCount;
   int lda = m, ldb = k, ldc = m;
   const real alpha = 1;
   const real beta = 0;

   blas::multMatrices(CBLAS_ORDER::CblasColMajor,
               CBLAS_TRANSPOSE::CblasNoTrans,
               CBLAS_TRANSPOSE::CblasNoTrans,
               m, n, k,
               alpha,
               b, ldb, a, lda,
               beta,
               result, ldc);
}