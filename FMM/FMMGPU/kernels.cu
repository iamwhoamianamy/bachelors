#include "kernels.cuh"
#include <stdio.h>
#include <thrust/complex.h>

namespace kernels
{
   __global__ void translateAllGPUKernelSimple(
      Vector3* result, const real* a, const Vector3* b,
      size_t harmonicCount, size_t harmonicOrder)
   {
      size_t harmonicLength = (harmonicOrder + 1) * (harmonicOrder + 1);
      uint harmonicId = blockIdx.x * blockDim.x + threadIdx.x;

      if(harmonicId < harmonicCount)
      {
         size_t harmonicBegin = harmonicLength * harmonicId;

         for(int l = 0; l < harmonicOrder; l++)
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
      size_t harmonicCount, size_t harmonicOrder)
   {
      size_t harmonicLength = (harmonicOrder + 1) * (harmonicOrder + 1);
      uint harmonicId = blockIdx.x * blockDim.x + threadIdx.x;

      if(harmonicId < harmonicCount)
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
                  zeroRes += b[lmToIndex(harmonicBegin, lambda, mu)] *
                     a[lmToIndex(harmonicBegin, dl, mu)] *
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
      size_t harmonicCount, size_t harmonicOrder)
   {
      size_t harmonicBegin = (harmonicOrder + 1) * (harmonicOrder + 1) * blockIdx.x;
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
      size_t harmonicCount, size_t harmonicOrder)
   {
      __shared__ real localA[kernels::HARMONIC_LENGTH];
      __shared__ Vector3 localB[kernels::HARMONIC_LENGTH];

      size_t harmonicBegin = blockIdx.x * (harmonicOrder + 1) * (harmonicOrder + 1);

      localA[threadIdx.x] = a[harmonicBegin + threadIdx.x];
      localB[threadIdx.x] = b[harmonicBegin + threadIdx.x];
      __syncthreads();

      if(threadIdx.x < harmonicOrder)
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

            result[lmToIndex(harmonicBegin, l, m)] = realRes * math::R_SQRT_2;
            result[lmToIndex(harmonicBegin, l, -m)] = imagRes * math::R_SQRT_2;
         }
      }
   }

   __global__ void matMulKernel(
      ComplexKernelMatrix A,
      ComplexKernelMatrix B,
      ComplexKernelMatrix C)
   {
      // Block row and column
      int blockRow = blockIdx.y;
      int blockCol = blockIdx.x;

      // Each thread block computes one sub-matrix Csub of C
      ComplexKernelMatrix Csub = getSubMatrix(C, blockRow, blockCol);

      // Each thread computes one element of Csub
      // by accumulating results into Cvalue
      thrust::complex<real> Cvalue = 0.0;

      // Thread row and column within Csub
      int row = threadIdx.y;
      int col = threadIdx.x;

      // Loop over all the sub-matrices of A and B that are
      // required to compute Csub
      // Multiply each pair of sub-matrices together
      // and accumulate the results
      for(int m = 0; m < (A.width / THREADS_PER_BLOCK); ++m)
      {

         // Get sub-matrix Asub of A
         ComplexKernelMatrix Asub = getSubMatrix(A, blockRow, m);

         // Get sub-matrix Bsub of B
         ComplexKernelMatrix Bsub = getSubMatrix(B, m, blockCol);

         // Shared memory used to store Asub and Bsub respectively
         __shared__ thrust::complex<real> As[THREADS_PER_BLOCK][THREADS_PER_BLOCK];
         __shared__ thrust::complex<real> Bs[THREADS_PER_BLOCK][THREADS_PER_BLOCK];

         // Load Asub and Bsub from device memory to shared memory
         // Each thread loads one element of each sub-matrix
         As[row][col] = getElement(Asub, row, col);
         Bs[row][col] = getElement(Bsub, row, col);

         // Synchronize to make sure the sub-matrices are loaded
         // before starting the computation
         __syncthreads();
         // Multiply Asub and Bsub together
         for(int e = 0; e < THREADS_PER_BLOCK; ++e)
            Cvalue = Cvalue + As[row][e] * Bs[e][col];

         // Synchronize to make sure that the preceding
         // computation is done before loading two new
         // sub-matrices of A and B in the next iteration
         __syncthreads();
      }

      // Write Csub to device memory
      // Each thread writes one element
      setElement(Csub, row, col, Cvalue);
   }
}