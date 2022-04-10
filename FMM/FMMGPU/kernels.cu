#include "kernels.cuh"
#include <stdio.h>

namespace kernels
{
   __global__ void translateAllGPUKernelSimple(
      Vector3* result, const real* a, const Vector3* b,
      size_t count, size_t order)
   {
      size_t harmonicLength = (order + 1) * (order + 1);
      uint harmonicId = blockIdx.x * blockDim.x + threadIdx.x;

      if(harmonicId < count)
      {
         size_t harmonicBegin = harmonicLength * harmonicId;

         for(int l = 0; l < order; l++)
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
                        strangeFactor(0, mu);
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

                        realRes += (RR * RM - IR * IM) * strangeFactor(m, mu);
                        imagRes += (RR * IM + IR * RM) * strangeFactor(m, mu);
                     }

                     if(-dl <= dnm && dnm <= dl)
                     {
                        RnM = getReal(a, harmonicBegin, dl, dnm);
                        InM = getImag(a, harmonicBegin, dl, dnm);

                        realRes += (RR * RnM - IR * InM) * strangeFactor(-m, mu);
                        imagRes -= (RR * InM + IR * RnM) * strangeFactor(-m, mu);
                     }
                  }
               }

               result[lmToIndex(harmonicBegin, l, m)] = realRes * math::R_SQRT_2;
               result[lmToIndex(harmonicBegin, l, -m)] = imagRes * math::R_SQRT_2;
            }
         }
      }
   }

   __global__ void translateAllGPUKernelSimpleXY(
      Vector3* result, const real* a, const Vector3* b,
      size_t count, size_t order)
   {
      size_t harmonicLength = (order + 1) * (order + 1);
      uint harmonicId = blockIdx.x * blockDim.x + threadIdx.x;

      if(harmonicId < count)
      {
         size_t harmonicBegin = harmonicLength * harmonicId;

         Vector3 zeroRes = 0;
         size_t l = threadIdx.y;

         for(int lambda = 0; lambda <= l; lambda++)
         {
            int dl = l - lambda;

            for(int mu = -lambda; mu <= lambda; mu++)
            {
               if(-dl <= mu && mu <= +dl)
               {
                  zeroRes += getHarmonic(b, harmonicBegin, lambda, mu) *
                     getHarmonic(a, harmonicBegin, dl, mu) *
                     strangeFactor(0, mu);
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

                     realRes += (RR * RM - IR * IM) * strangeFactor(m, mu);
                     imagRes += (RR * IM + IR * RM) * strangeFactor(m, mu);
                  }

                  if(-dl <= dnm && dnm <= dl)
                  {
                     RnM = getReal(a, harmonicBegin, dl, dnm);
                     InM = getImag(a, harmonicBegin, dl, dnm);

                     realRes += (RR * RnM - IR * InM) * strangeFactor(-m, mu);
                     imagRes -= (RR * InM + IR * RnM) * strangeFactor(-m, mu);
                  }
               }
            }

            result[lmToIndex(harmonicBegin, l, m)] = realRes * math::R_SQRT_2;
            result[lmToIndex(harmonicBegin, l, -m)] = imagRes * math::R_SQRT_2;
         }
      }
   }

   __global__ void translateAllGPUKernelBlockForHarmonic(
      Vector3* result, const real* a, const Vector3* b,
      size_t count, size_t order)
   {
      size_t harmonicBegin = (order + 1) * (order + 1) * blockIdx.x;
      Vector3 zeroRes = 0;
      int l = threadIdx.x;

      for(int lambda = 0; lambda <= l; lambda++)
      {
         int dl = l - lambda;

         for(int mu = -lambda; mu <= lambda; mu++)
         {
            if(-dl <= mu && mu <= +dl)
            {
               zeroRes += getHarmonic(b, harmonicBegin, lambda, mu) *
                  getHarmonic(a, harmonicBegin, dl, mu) *
                  strangeFactor(0, mu);
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

                  realRes += (RR * RM - IR * IM) * strangeFactor(m, mu);
                  imagRes += (RR * IM + IR * RM) * strangeFactor(m, mu);
               }

               if(-dl <= dnm && dnm <= dl)
               {
                  RnM = getReal(a, harmonicBegin, dl, dnm);
                  InM = getImag(a, harmonicBegin, dl, dnm);

                  realRes += (RR * RnM - IR * InM) * strangeFactor(-m, mu);
                  imagRes -= (RR * InM + IR * RnM) * strangeFactor(-m, mu);
               }
            }
         }

         result[lmToIndex(harmonicBegin, l, m)] = realRes * math::R_SQRT_2;
         result[lmToIndex(harmonicBegin, l, -m)] = imagRes * math::R_SQRT_2;
      }
   }

   __global__ void translateAllGPUKernelBlockForHarmonicShared(
      Vector3* result, const real* a, const Vector3* b,
      size_t count, size_t order)
   {
      __shared__ real localA[kernels::HARMONIC_LENGTH];
      __shared__ Vector3 localB[kernels::HARMONIC_LENGTH];
      __shared__ Vector3 localRes[kernels::HARMONIC_LENGTH];

      size_t harmonicBegin = blockIdx.x * (order + 1) * (order + 1);

      localA[threadIdx.x] = a[harmonicBegin + threadIdx.x];
      __syncthreads();

      localB[threadIdx.x] = b[harmonicBegin + threadIdx.x];
      __syncthreads();

      localRes[threadIdx.x] = 0;
      __syncthreads();

      if(threadIdx.x < order)
      {
         int l = threadIdx.x;
         Vector3 zeroRes = 0;

         for(int lambda = 0; lambda <= l; lambda++)
         {
            int dl = l - lambda;

            for(int mu = -lambda; mu <= lambda; mu++)
            {
               if(-dl <= mu && mu <= +dl)
               {
                  zeroRes += getHarmonic(localB, lambda, mu) *
                     getHarmonic(localA, dl, mu) *
                     strangeFactor(0, mu);
               }
            }
         }

         localRes[lmToIndex(l, 0)] = zeroRes;

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

                  Vector3 RR = getReal(localB, lambda, mu);
                  Vector3 IR = getImag(localB, lambda, mu);

                  real RM = 0;
                  real IM = 0;

                  real RnM = 0;
                  real InM = 0;

                  if(-dl <= dm && dm <= dl)
                  {
                     RM = getReal(localA, dl, dm);
                     IM = getImag(localA, dl, dm);

                     realRes += (RR * RM - IR * IM) * strangeFactor(m, mu);
                     imagRes += (RR * IM + IR * RM) * strangeFactor(m, mu);
                  }

                  if(-dl <= dnm && dnm <= dl)
                  {
                     RnM = getReal(localA, dl, dnm);
                     InM = getImag(localA, dl, dnm);

                     realRes += (RR * RnM - IR * InM) * strangeFactor(-m, mu);
                     imagRes -= (RR * InM + IR * RnM) * strangeFactor(-m, mu);
                  }
               }
            }

            localRes[lmToIndex(l, m)] = realRes * math::R_SQRT_2;
            localRes[lmToIndex(l, -m)] = imagRes * math::R_SQRT_2;
         }
      }

      __syncthreads();

      result[harmonicBegin + threadIdx.x] = localRes[threadIdx.x];
   }
}