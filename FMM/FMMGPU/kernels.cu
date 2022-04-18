#include "kernels.cuh"
#include <stdio.h>
#include <thrust/complex.h>

namespace kernels
{
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

   __global__ void matMulKernel(
      ComplexKernelMatrix A,
      ComplexKernelMatrix B,
      ComplexKernelMatrix C)
   {
   //   int blockRow = blockIdx.y;
   //   int blockCol = blockIdx.x;

   //   ComplexKernelMatrix Csub = getSubMatrix(C, blockRow, blockCol);
   //   Complex Cvalue = make_cuComplex(0, 0);

   //   int row = threadIdx.y;
   //   int col = threadIdx.x;

   //   for(int m = 0; m < (A.width / THREADS_PER_BLOCK); ++m)
   //   {
   //      ComplexKernelMatrix Asub = getSubMatrix(A, blockRow, m);
   //      ComplexKernelMatrix Bsub = getSubMatrix(B, m, blockCol);

   //      __shared__ Complex As[THREADS_PER_BLOCK][THREADS_PER_BLOCK];
   //      __shared__ Complex Bs[THREADS_PER_BLOCK][THREADS_PER_BLOCK];

   //      As[row][col] = getElement(Asub, row, col);
   //      Bs[row][col] = getElement(Bsub, row, col);

   //      __syncthreads();

   //      for(int e = 0; e < THREADS_PER_BLOCK; ++e)
   //         Cvalue = Cvalue + As[row][e] * Bs[e][col];

   //      __syncthreads();
   //   }

   //   setElement(Csub, row, col, Cvalue);
   }
}