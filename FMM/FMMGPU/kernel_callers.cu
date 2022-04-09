#pragma once
#include "kernel_callers.hpp"
#include <vector>
#include "kernels.cuh"
#include "cuda_helper.hpp"
#include "dev_ptr.hpp"

namespace kernels
{
   std::vector<real> addVectors(const std::vector<real>& a,
                                const std::vector<real>& b)
   {
      size_t size = a.size();

      cuda::DevPtr<real> dev_a(a.data(), a.size());
      cuda::DevPtr<real> dev_b(b.data(), b.size());
      cuda::DevPtr<real> dev_res(size);

      addingKernel<<<1, size>>>(dev_res.data(), dev_a.data(), dev_b.data());

      cuda::tryKernelLaunch();
      cuda::tryKernelSynchronize();

      real* res = new real[size];
      dev_res.copyToHost(res);

      return std::vector<real>(res, res + size);
   }

   void translateAllCPU(Vector3* result,
                        const real* a,
                        const Vector3* b,
                        size_t count, size_t order)
   {
      size_t harmonicLength = (order + 1) * (order + 1);

      for(size_t harmonicId = 0; harmonicId < count; harmonicId++)
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
                     zeroRes.x += getHarmonic(b, harmonicBegin, lambda, mu)->x *
                        *getHarmonic(a, harmonicBegin, dl, mu) *
                        Harmonics::strangeFactor(0, mu);
                     zeroRes.y += getHarmonic(b, harmonicBegin, lambda, mu)->y *
                        *getHarmonic(a, harmonicBegin, dl, mu) *
                        Harmonics::strangeFactor(0, mu);
                     zeroRes.z += getHarmonic(b, harmonicBegin, lambda, mu)->z *
                        *getHarmonic(a, harmonicBegin, dl, mu) *
                        Harmonics::strangeFactor(0, mu);
                  }
               }
            }

            getHarmonic(result, harmonicBegin, l, 0)->x = zeroRes.x;
            getHarmonic(result, harmonicBegin, l, 0)->y = zeroRes.y;
            getHarmonic(result, harmonicBegin, l, 0)->z = zeroRes.z;

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

                        realRes.x += (RR.x * RM - IR.x * IM) * Harmonics::strangeFactor(m, mu);
                        realRes.y += (RR.y * RM - IR.y * IM) * Harmonics::strangeFactor(m, mu);
                        realRes.z += (RR.z * RM - IR.z * IM) * Harmonics::strangeFactor(m, mu);

                        imagRes.x += (RR.x * IM + IR.x * RM) * Harmonics::strangeFactor(m, mu);
                        imagRes.y += (RR.y * IM + IR.y * RM) * Harmonics::strangeFactor(m, mu);
                        imagRes.z += (RR.z * IM + IR.z * RM) * Harmonics::strangeFactor(m, mu);
                     }

                     if(-dl <= dnm && dnm <= dl)
                     {
                        RnM = getReal(a, harmonicBegin, dl, dnm);
                        InM = getImag(a, harmonicBegin, dl, dnm);

                        realRes.x += (RR.x * RnM - IR.x * InM) * Harmonics::strangeFactor(-m, mu);
                        realRes.y += (RR.y * RnM - IR.y * InM) * Harmonics::strangeFactor(-m, mu);
                        realRes.z += (RR.z * RnM - IR.z * InM) * Harmonics::strangeFactor(-m, mu);

                        imagRes.x -= (RR.x * InM + IR.x * RnM) * Harmonics::strangeFactor(-m, mu);
                        imagRes.y -= (RR.y * InM + IR.y * RnM) * Harmonics::strangeFactor(-m, mu);
                        imagRes.z -= (RR.z * InM + IR.z * RnM) * Harmonics::strangeFactor(-m, mu);
                     }
                  }
               }

               getHarmonic(result, harmonicBegin, l, m)->x = realRes.x * math::R_SQRT_2;

               getHarmonic(result, harmonicBegin, l, -m)->x = imagRes.x * math::R_SQRT_2;
               getHarmonic(result, harmonicBegin, l, -m)->y = imagRes.y * math::R_SQRT_2;
               getHarmonic(result, harmonicBegin, l, -m)->z = imagRes.z * math::R_SQRT_2;
            }
         }
      }
   }
}