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

   void translateAllGPU(Vector3* result,
                        const real* a,
                        const Vector3* b,
                        size_t count, size_t order)
   {
      uint BLOCK_COUNT = (count + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
      size_t harmonicLength = (order + 1) * (order + 1);

      cuda::DevPtr<Vector3> result_dev(count * harmonicLength);
      cuda::DevPtr<real> a_dev(a, count * harmonicLength, 0);
      cuda::DevPtr<Vector3> b_dev(b, count * harmonicLength, 0);

      kernels::translateAllGPUKernel<<<BLOCK_COUNT, THREADS_PER_BLOCK>>>
         (result_dev.data(), a_dev.data(), b_dev.data(), count, order);

      cuda::tryKernelLaunch();
      cuda::tryKernelSynchronize();

      result_dev.copyToHost(result);
   }

   size_t lmToIndex(int harmonicBegin,
                    int l, int m)
   {
      return harmonicBegin + l * l + l + m;
   }

   __all__ real strangeFactor(int m, int mu)
   {
      return pow(-1, -0.5 * (abs(m) - abs(mu) - abs(m - mu)));
   }
}