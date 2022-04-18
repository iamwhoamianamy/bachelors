#pragma once
#include "translation_algorithms.hpp"
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

   void translateAllGPU(Vector3* result,
                        const real* a,
                        const Vector3* b,
                        size_t harmonicCount, size_t harmonicOrder)
   {
      size_t harmonicLength = (harmonicOrder + 1) * (harmonicOrder + 1);

      cuda::DevPtr<Vector3> result_dev(harmonicCount * harmonicLength);
      cuda::DevPtr<real> a_dev(a, harmonicCount * harmonicLength, 0);
      cuda::DevPtr<Vector3> b_dev(b, harmonicCount * harmonicLength, 0);

      dim3 BLOCKS((harmonicCount + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
      dim3 THREADS(THREADS_PER_BLOCK, harmonicOrder);

      kernels::translateAllGPUKernelSimpleXY<<<BLOCKS, THREADS>>>
         (result_dev.data(), a_dev.data(), b_dev.data(), harmonicCount, harmonicOrder);
      
      cuda::tryKernelLaunch();
      cuda::tryKernelSynchronize();

      result_dev.copyToHost(result);
   }

   void kernels::translateAllGPUMatrix(
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

      cuda::DevPtr<Complex> aDev(a, harCountPadded * harLenPadded);
      ComplexKernelMatrix A;

      A.width = A.stride = harLenPadded;
      A.height = harCountPadded;
      A.elements = aDev.data();

      cuda::DevPtr<Complex> bDev(b, harLenPadded * harLenPadded);
      ComplexKernelMatrix B;

      B.width = B.stride = harLenPadded;
      B.height = harLenPadded;
      B.elements = bDev.data();

      cuda::DevPtr<Complex> cDev(harCountPadded * harLenPadded);
      ComplexKernelMatrix C;

      C.width = C.stride = harLenPadded;
      C.height = harCountPadded;
      C.elements = cDev.data();

      dim3 dimBlock(THREADS_PER_BLOCK, THREADS_PER_BLOCK);
      dim3 dimGrid(B.width / dimBlock.x, A.height / dimBlock.y);
      matMulKernel<<<dimGrid, dimBlock >>>(A, B, C);

      cuda::tryKernelLaunch();
      cuda::tryKernelSynchronize();

      cDev.copyToHost(result);
   }

   size_t lmToIndex(int harmonicBegin,
                    int l, int m)
   {
      return harmonicBegin + l * l + l + m;
   }

   size_t lmToIndex(int l, int m)
   {
      return l * l + l + m;
   }

   __all__ real strangeFactor(int m, int mu)
   {
      return pow(-1, -0.5 * (abs(m) - abs(mu) - abs(m - mu)));
   }
}