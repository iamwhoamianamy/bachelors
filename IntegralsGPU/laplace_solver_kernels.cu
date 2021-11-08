#include "laplace_solver_kernels.cuh";
#include "laplace_data.cuh"

typedef unsigned int uint;

#if (!defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600)
#else
__device__ double atomicAdd(double* a, double b) { return b; }
#endif

__global__ void laplace_solver_kernels::SolverKernelArrays(
   const double* quadratures_X,
   const double* quadratures_Y,
   const double* quadratures_Z,
   const double* normals_X,
   const double* normals_Y,
   const double* normals_Z,
   const double* points_X,
   const double* points_Y,
   const double* points_Z,
   const double* weights,
   const double* areas,
   const int trianglesCount,
   const int pointsCount,
   const int quadraturesOrder,
   double* result)
{
   uint p = blockIdx.x;
   extern __shared__ double shared_array[];
   double* integral = (__shared__ double*)shared_array;
   uint t = threadIdx.x;

   while(t < trianglesCount)
   {
      double tringle_sum_1 = 0;
      double tringle_sum_2 = 0;

      for(size_t o = 0; o < quadraturesOrder; o++)
      {
         int ind = t * quadraturesOrder + o;
         tringle_sum_1 += weights[o] * laplace_data::laplaceIntegral1(
            quadratures_X[ind], quadratures_Y[ind], quadratures_Z[ind],
            points_X[p], points_Y[p], points_Z[p],
            normals_X[t], normals_Y[t], normals_Z[t]);

         tringle_sum_2 += weights[o] * laplace_data::laplaceIntegral2(
            quadratures_X[ind], quadratures_Y[ind], quadratures_Z[ind],
            points_X[p], points_Y[p], points_Z[p],
            normals_X[t], normals_Y[t], normals_Z[t]);
      }

      atomicAdd(&integral[t % blockDim.x], (tringle_sum_1 - tringle_sum_2) * areas[t]);
      t += blockDim.x;
   }

   __syncthreads();
   t = threadIdx.x;

   for(unsigned int s = 1; s < blockDim.x; s *= 2)
   {
      if(t % (2 * s) == 0)
      {
         if(t + s < trianglesCount)
         {
            integral[t] += integral[t + s];
         }
      }

      __syncthreads();
   }

   __syncthreads();
   t = threadIdx.x;

   if(t == 0)
   {
      result[p] = integral[0] / (4.0 * laplace_data::PI);
   }
}

__global__ void laplace_solver_kernels::SolverKernelVector3s(
   const Vector3* quadratures,
   const Vector3* normals,
   const Vector3* points,
   const double* weights,
   const double* areas,
   const int trianglesCount,
   const int pointsCount,
   const int quadraturesOrder,
   double* result)
{
   uint p = blockIdx.x;
   extern __shared__ double shared_array[];
   double* integral = (__shared__ double*)shared_array;
   uint t = threadIdx.x;

   while(t < trianglesCount)
   {
      double tringle_sum_1 = 0;
      double tringle_sum_2 = 0;

      for(size_t o = 0; o < quadraturesOrder; o++)
      {
         int ind = t * quadraturesOrder + o;
         tringle_sum_1 += weights[o] * laplace_data::laplaceIntegral1(
            quadratures[ind], points[p], normals[t]);

         tringle_sum_2 += weights[o] * laplace_data::laplaceIntegral2(
            quadratures[ind], points[p], normals[t]);
      }

      atomicAdd(&integral[t % blockDim.x], (tringle_sum_1 - tringle_sum_2) * areas[t]);
      t += blockDim.x;
   }

   __syncthreads();
   t = threadIdx.x;

   for(unsigned int s = 1; s < blockDim.x; s *= 2)
   {
      if(t % (2 * s) == 0)
      {
         if(t + s < trianglesCount)
         {
            integral[t] += integral[t + s];
         }
      }

      __syncthreads();
   }

   __syncthreads();
   t = threadIdx.x;

   if(t == 0)
   {
      result[p] = integral[0] / (4.0 * laplace_data::PI);
   }
}