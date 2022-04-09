#include "kernels.cuh"
#include <stdio.h>

namespace kernels
{
   __global__ void translateAllGPUKernel(Vector3* result,
                              const real* a,
                              const Vector3* b,
                              size_t count, size_t order)
   {
      size_t harmonicLength = (order + 1) * (order + 1);
      
      uint harmonicId = blockIdx.x * blockDim.x + threadIdx.x;

      if(harmonicId < count)
      {
         size_t harmonicBegin = harmonicLength * harmonicId;
         size_t harmonicEnd = harmonicBegin + harmonicLength;

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
}