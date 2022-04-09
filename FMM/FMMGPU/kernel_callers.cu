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
            real zeroResX = 0;
            real zeroResY = 0;
            real zeroResZ = 0;

            for(int lambda = 0; lambda <= l; lambda++)
            {
               int dl = l - lambda;

               for(int mu = -lambda; mu <= lambda; mu++)
               {
                  if(-dl <= mu && mu <= +dl)
                  {
                     zeroResX += getHarmonic(b, harmonicBegin, lambda, mu)->x *
                        *getHarmonic(a, harmonicBegin, dl, mu) *
                        Harmonics::strangeFactor(0, mu);
                     zeroResY += getHarmonic(b, harmonicBegin, lambda, mu)->y *
                        *getHarmonic(a, harmonicBegin, dl, mu) *
                        Harmonics::strangeFactor(0, mu);
                     zeroResZ += getHarmonic(b, harmonicBegin, lambda, mu)->z *
                        *getHarmonic(a, harmonicBegin, dl, mu) *
                        Harmonics::strangeFactor(0, mu);
                  }
               }
            }

            getHarmonic(result, harmonicBegin, l, 0)->x =zeroResX;
            getHarmonic(result, harmonicBegin, l, 0)->y =zeroResY;
            getHarmonic(result, harmonicBegin, l, 0)->z =zeroResZ;

            for(int m = 1; m <= l; m++)
            {
               real realResX = 0;
               real realResY = 0;
               real realResZ = 0;

               real imagResX = 0;
               real imagResY = 0;
               real imagResZ = 0;

               for(int lambda = 0; lambda <= l; lambda++)
               {
                  int dl = l - lambda;

                  for(int mu = -lambda; mu <= lambda; mu++)
                  {
                     int dm = m - mu;
                     int dnm = -m - mu;

                     real RRx = getReal(b, harmonicBegin, lambda, mu).x;
                     real RRy = getReal(b, harmonicBegin, lambda, mu).y;
                     real RRz = getReal(b, harmonicBegin, lambda, mu).z;

                     real IRx = getImag(b, harmonicBegin, lambda, mu).x;
                     real IRy = getImag(b, harmonicBegin, lambda, mu).y;
                     real IRz = getImag(b, harmonicBegin, lambda, mu).z;

                     real RM = 0;
                     real IM = 0;

                     real RnM = 0;
                     real InM = 0;

                     if(-dl <= dm && dm <= dl)
                     {
                        RM = getReal(a, harmonicBegin, dl, dm);
                        IM = getImag(a, harmonicBegin, dl, dm);

                        realResX += (RRx * RM - IRx * IM) * Harmonics::strangeFactor(m, mu);
                        realResY += (RRy * RM - IRy * IM) * Harmonics::strangeFactor(m, mu);
                        realResZ += (RRz * RM - IRz * IM) * Harmonics::strangeFactor(m, mu);

                        imagResX += (RRx * IM + IRx * RM) * Harmonics::strangeFactor(m, mu);
                        imagResY += (RRy * IM + IRy * RM) * Harmonics::strangeFactor(m, mu);
                        imagResZ += (RRz * IM + IRz * RM) * Harmonics::strangeFactor(m, mu);
                     }

                     if(-dl <= dnm && dnm <= dl)
                     {
                        RnM = getReal(a, harmonicBegin, dl, dnm);
                        InM = getImag(a, harmonicBegin, dl, dnm);

                        realResX += (RRx * RnM - IRx * InM) * Harmonics::strangeFactor(-m, mu);
                        realResY += (RRy * RnM - IRy * InM) * Harmonics::strangeFactor(-m, mu);
                        realResZ += (RRz * RnM - IRz * InM) * Harmonics::strangeFactor(-m, mu);

                        imagResX -= (RRx * InM + IRx * RnM) * Harmonics::strangeFactor(-m, mu);
                        imagResY -= (RRy * InM + IRy * RnM) * Harmonics::strangeFactor(-m, mu);
                        imagResZ -= (RRz * InM + IRz * RnM) * Harmonics::strangeFactor(-m, mu);
                     }
                  }
               }

               getHarmonic(result, harmonicBegin, l, m)->x = realResX * math::R_SQRT_2;
               getHarmonic(result, harmonicBegin, l, m)->y = realResY * math::R_SQRT_2;
               getHarmonic(result, harmonicBegin, l, m)->z = realResZ * math::R_SQRT_2;

               getHarmonic(result, harmonicBegin, l, -m)->x = imagResX * math::R_SQRT_2;
               getHarmonic(result, harmonicBegin, l, -m)->y = imagResY * math::R_SQRT_2;
               getHarmonic(result, harmonicBegin, l, -m)->z = imagResZ * math::R_SQRT_2;
            }
         }
      }
   }
}