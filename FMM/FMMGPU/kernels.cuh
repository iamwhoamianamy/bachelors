#include "cuda_runtime.h"
#include "cublas_v2.h"
#include "device_launch_parameters.h"
#include "harmonics.hpp"

namespace kernels
{
   template<class T>
   __global__ void addingKernel(T* c, const T* a, const T* b)
   {
      uint i = threadIdx.x;
      c[i] = a[i] + b[i];
   }

   __global__ void translateAllGPUKernelSimple(
      Vector3* result, const real* a, const Vector3* b,
      size_t harmonicCount, size_t harmonicOrder);

   __global__ void translateAllGPUKernelSimpleXY(
      Vector3* result, const real* a, const Vector3* b,
      size_t harmonicCount, size_t harmonicOrder);

   __global__ void translateAllGPUKernelBlockForHarmonic(
      Vector3* result, const real* a, const Vector3* b,
      size_t harmonicCount, size_t harmonicOrder);

   __global__ void translateAllGPUKernelBlockForHarmonicShared(
      Vector3* result, const real* a, const Vector3* b,
      size_t harmonicCount, size_t harmonicOrder);

   template <class T>
   struct KernelMatrix
   {
      int width;
      int height;
      int stride;
      T* elements;
   };

   typedef KernelMatrix<thrust::complex<real>> ComplexKernelMatrix;

   template<class T>
   __device__ T getElement(const KernelMatrix<T> A, int row, int col)
   {
      return A.elements[row * A.stride + col];
   }

   template<class T>
   __device__ void setElement(KernelMatrix<T> A, int row, int col,
                              T value)
   {
      A.elements[row * A.stride + col] = value;
   }

   template<class T>
   __device__ KernelMatrix<T> getSubMatrix(KernelMatrix<T> A, int row, int col)
   {
      KernelMatrix<T> Asub;
      Asub.width = THREADS_PER_BLOCK;
      Asub.height = THREADS_PER_BLOCK;
      Asub.stride = A.stride;
      Asub.elements = &A.elements[A.stride * THREADS_PER_BLOCK * row
         + THREADS_PER_BLOCK * col];
      return Asub;
   }

   __global__ void matMulKernel(
      ComplexKernelMatrix A,
      ComplexKernelMatrix B,
      ComplexKernelMatrix C);

   __all__ size_t lmToIndex(int harmonicBegin,
                    int l, int m);

   __all__ size_t lmToIndex(int l, int m);

   template <typename T>
   inline __all__ T& getHarmonic(T* harmonics,
                  int harmonicBegin,
                  int l, int m)
   {
      return harmonics[harmonicBegin + l * l + l + m];
   }

   template <typename T>
   inline __all__ T& getHarmonic(T* harmonics,
                          int l, int m)
   {
      return harmonics[l * l + l + m];
   }

   template <typename T>
   inline __all__ T getReal(const T* harmonics,
             int harmonicBegin,
             int l, int m)
   {
      return harmonics[lmToIndex(harmonicBegin, l, abs(m))] *
         (math::R_SQRT_2 * (m != 0) + (m == 0));
   }

   template <typename T>
   inline __all__ T getReal(const T* harmonics,
                     int l, int m)
   {
      return harmonics[lmToIndex(l, abs(m))] *
         (math::R_SQRT_2 * (m != 0) + (m == 0));
   }

   template <typename T>
   inline __all__ T getImag(const T* harmonics,
             int harmonicBegin,
             int l, int m)
   {
     return harmonics[lmToIndex(harmonicBegin, l, -abs(m))] *
         math::R_SQRT_2 * (m != 0) * math::sign(m);
   }

   template <typename T>
   inline __all__ T getImag(const T* harmonics,
                     int l, int m)
   {
      return harmonics[lmToIndex(l, -abs(m))] *
         math::R_SQRT_2 * (m != 0) * math::sign(m);
   }

   __all__ real strangeFactor(int m, int mu);

   const uint THREADS_PER_BLOCK = 32;
   const uint HARMONIC_ORDER = 10;
   const uint HARMONIC_LENGTH = 121;
}