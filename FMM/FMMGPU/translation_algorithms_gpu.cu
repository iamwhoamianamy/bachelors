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

      {
         //dim3 BLOCKS((harmonicCount + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
         //dim3 THREADS(THREADS_PER_BLOCK);

         //kernels::translateAllGPUKernelSimple<<<BLOCKS, THREADS>>>
         //   (result_dev.data(), a_dev.data(), b_dev.data(), harmonicCount, harmonicOrder);
      }

      {
         dim3 BLOCKS((harmonicCount + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
         dim3 THREADS(THREADS_PER_BLOCK, harmonicOrder);

         kernels::translateAllGPUKernelSimpleXY<<<BLOCKS, THREADS>>>
            (result_dev.data(), a_dev.data(), b_dev.data(), harmonicCount, harmonicOrder);
      }

      {
      //   dim3 BLOCKS(harmonicCount);
      //   dim3 THREADS(harmonicOrder);

      //   kernels::translateAllGPUKernelBlockForHarmonic<<<BLOCKS, THREADS>>>
      //      (result_dev.data(), a_dev.data(), b_dev.data(), harmonicCount, harmonicOrder);
      }

      {
         //dim3 BLOCKS(harmonicCount);
         //dim3 THREADS(harmonicLength);

         //kernels::translateAllGPUKernelBlockForHarmonicShared<<<BLOCKS, THREADS>>>
         //   (result_dev.data(), a_dev.data(), b_dev.data(), harmonicCount, harmonicOrder);
      }

      cuda::tryKernelLaunch();
      cuda::tryKernelSynchronize();

      result_dev.copyToHost(result);
   }

   void translateAllGPUMatrix(
      Complex* result,
      const Complex* a,
      const Complex* b,
      size_t harmonicCount,
      size_t harmonicOrder)
   {     
      size_t harmonicLength = (harmonicOrder + 1) * (harmonicOrder + 1);
      size_t harLenPadded = math::nextDevisible(harmonicLength, THREADS_PER_BLOCK);
      size_t harCountPadded = math::nextDevisible(harmonicCount, THREADS_PER_BLOCK);

      ComplexKernelMatrix d_A;
      d_A.width = d_A.stride = harLenPadded;
      d_A.height = harCountPadded;
      size_t size = harCountPadded * harLenPadded * sizeof(Complex);
      cudaMalloc(&d_A.elements, size);
      cudaMemcpy(d_A.elements, a, size,
                 cudaMemcpyHostToDevice);

      ComplexKernelMatrix d_B;
      d_B.width = d_B.stride = harLenPadded;
      d_B.height = harLenPadded;
      size = harLenPadded * harLenPadded * sizeof(Complex);
      cudaMalloc(&d_B.elements, size);
      cudaMemcpy(d_B.elements, b, size,
                 cudaMemcpyHostToDevice);

      ComplexKernelMatrix d_C;
      d_C.width = d_C.stride = harLenPadded;
      d_C.height = harCountPadded;
      size = harCountPadded * harLenPadded * sizeof(Complex);
      cudaMalloc(&d_C.elements, size);

      // Invoke kernel
      dim3 dimBlock(THREADS_PER_BLOCK, THREADS_PER_BLOCK);
      dim3 dimGrid(d_B.width / dimBlock.x, d_A.height / dimBlock.y);
      matMulKernel<<<dimGrid, dimBlock >>>(d_A, d_B, d_C);

      cuda::tryKernelLaunch();
      cuda::tryKernelSynchronize();

      // Read C from device memory
      cudaMemcpy(result, d_C.elements, size,
                 cudaMemcpyDeviceToHost);

      // Free device memory
      cudaFree(d_A.elements);
      cudaFree(d_B.elements);
      cudaFree(d_C.elements);
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