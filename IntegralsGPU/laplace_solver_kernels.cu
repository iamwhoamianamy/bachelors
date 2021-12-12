#include "laplace_solver_kernels.cuh";
#include "laplace_data.cuh"

typedef unsigned int uint;

#if (!defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600)
#else
__device__ real atomicAdd(double* a, double b) { return b; }
#endif

using namespace laplace_data;

__global__ void laplace_solver_kernels::solverKernelVector3s(
   const Vector3* quadPoints, const Vector3* normals, const Vector3* points,
   const real* weights, const real* areas,
   const int trianglesCount, const int pointsCount, const int quadraturesOrder,
   real* results)
{
   // Slice of points in global memory stored in shared memory
   __shared__ Vector3 points_shared[POINTS_PER_BLOCK];

   // Index of point in global and shared memory
   uint point_idx_g = threadIdx.x + POINTS_PER_BLOCK * blockIdx.x;
   uint point_idx_l = threadIdx.x;

   // Loading points into shared memory
   points_shared[point_idx_l] = points[point_idx_g];

   // Result for each point(thread)
   __shared__ real results_shared[POINTS_PER_BLOCK];

   // Initializing results in shared memory
   results_shared[point_idx_l] = 0;

   __syncthreads();

   const int quadratures_count = trianglesCount * quadraturesOrder;

   // Iterating through every block of quadPoints
   for(size_t b = 0; b < quadratures_count / QUADS_PER_BLOCK; b++)
   {
      // Slice of quadPoints in global memory stored in shared memory 
      __shared__ Vector3 quads_shared[QUADS_PER_BLOCK];

      // Weights for current block
      __shared__ real weights_shared[QUADS_PER_BLOCK];

      // Normals for current block
      __shared__ Vector3 normals_shared[TR_PER_BLOCK];

      // Areas for current block
      __shared__ real areas_shared[TR_PER_BLOCK];

      // Loading normals and areas for triangles in current block
      if(threadIdx.x == 0)
      {
         for(size_t i = 0; i < TR_PER_BLOCK; i++)
         {
            normals_shared[i] = normals[b * TR_PER_BLOCK + i];
            areas_shared[i] = areas[b * TR_PER_BLOCK + i];
         }
      }

      // Index of quadrature to load from global to shared memory
      const uint quad_idx_g = threadIdx.x + QUADS_PER_BLOCK * b;
      const uint quad_idx_l = threadIdx.x;

      // Loading quadPoints into shared memory
      quads_shared[quad_idx_l] = quadPoints[quad_idx_g];
      weights_shared[quad_idx_l] = weights[quad_idx_l % quadraturesOrder];

      __syncthreads();

      // Iterating through every triangle
      for(size_t t = 0; t < TR_PER_BLOCK; t++)
      {
         real triangle_sum = 0.0;
         real triangle_sum_2 = 0.0;

         // Iterating through every quadrature in triangle
         for(size_t q = 0; q < quadraturesOrder; q++)
         {
            const int idx = t * quadraturesOrder + q;

            triangle_sum += weights_shared[idx] * calcLaplaceIntegralGPU(
               quads_shared[idx].x, quads_shared[idx].x, quads_shared[idx].x,
               points_shared[point_idx_l].x, points_shared[point_idx_l].y, points_shared[point_idx_l].z,
               normals_shared[t].x, normals_shared[t].y, normals_shared[t].z);
         }

         results_shared[point_idx_l] += (triangle_sum - triangle_sum_2) * areas_shared[t];
      }

      __syncthreads();
   }

   results[point_idx_g] = results_shared[point_idx_l] / (4.0 * PI);

#pragma region old
   //// Slice of points in global memory stored in shared memory
   //__shared__ Vector3 points_shared[POINTS_PER_BLOCK];

   //// Index of point in global and shared memory
   //const uint point_idx_g = threadIdx.y + blockDim.y * blockIdx.y;
   //const uint point_idx_l = threadIdx.y;

   //// Loading points into shared memory
   //if(threadIdx.x == 0)
   //   points_shared[point_idx_l] = points[point_idx_g];

   //// Result for each point(thread)
   //__shared__ real results_shared[POINTS_PER_BLOCK];

   //// Initializing results in shared memory
   //if(threadIdx.x == 0)
   //   results_shared[point_idx_l] = 0;

   //__syncthreads();

   //const int quadratures_count = trianglesCount * quadraturesOrder;
   //
   //// Iterating through every block of quadPoints
   //for(size_t b = 0; b < quadratures_count / QUADS_PER_BLOCK; b++)
   //{
   //   // Slice of quadPoints in global memory stored in shared memory 
   //   __shared__ Vector3 quads_shared[QUADS_PER_BLOCK];

   //   // Weights for current block
   //   __shared__ real weights_shared[QUADS_PER_BLOCK];

   //   // Normal for current block
   //   __shared__ Vector3 normals_shared;

   //   // Area for current block
   //   __shared__ real areas_shared;

   //   // Loading normals_shared for current block
   //   if(threadIdx.x == 0 && threadIdx.y == 0)
   //   {
   //      normals_shared = normals[b];
   //      areas_shared = areas[b];
   //   }

   //   // Index of quadrature to load from global to shared memory
   //   const uint quad_idx_g = threadIdx.x + blockDim.x * b;
   //   const uint quad_idx_l = threadIdx.x;

   //   // Loading quadPoints into shared memory
   //   if(threadIdx.y == 0)
   //   {
   //      quads_shared[quad_idx_l] = quadPoints[quad_idx_g];
   //      weights_shared[quad_idx_l] = weights[quad_idx_l];
   //   }
   //   /*if(quad_idx_g < quadratures_count)
   //   {
   //      
   //   }
   //   else
   //   {
   //      weights_shared[quad_idx_l] = 0;
   //   }*/

   //   __syncthreads();

   //   if(threadIdx.x == 0)
   //   {
   //      real block_sum_1 = 0.0;
   //      real block_sum_2 = 0.0;

   //      // Iterating through every quadrature in block
   //      for(size_t q = 0; q < QUADS_PER_BLOCK; q++)
   //      {
   //         block_sum_1 += weights_shared[q] * laplaceIntegral1GPU(
   //            quads_shared[q], points_shared[point_idx_l], normals_shared);

   //         block_sum_2 += weights_shared[q] * laplaceIntegral2GPU(
   //            quads_shared[q], points_shared[point_idx_l], normals_shared);
   //      }

   //      results_shared[point_idx_l] += (block_sum_1 - block_sum_2) * areas_shared;
   //   }

   //   __syncthreads();
   //}

   //if(threadIdx.x == 0)
   //   results[point_idx_g] = results_shared[point_idx_l] / (4.0 * PI);
#pragma endregion old
}

__global__ void laplace_solver_kernels::solverKernelStructs(
   const QuadPoint* quadPoints, const Vector3* points, const int quadraturesCount,
   real* results)
{
   __shared__ Vector3 points_shared[BLOCK_SIZE];
   __shared__ QuadPoint quads_shared[BLOCK_SIZE];

   const uint point_idx_g = threadIdx.x + BLOCK_SIZE * blockIdx.x;
   points_shared[threadIdx.x] = points[point_idx_g];

   real point_result = 0.0;

   for(size_t b = 0; b < (quadraturesCount + BLOCK_SIZE - 1) / BLOCK_SIZE; b++)
   {
      quads_shared[threadIdx.x] = quadPoints[threadIdx.x + BLOCK_SIZE * b];

      __syncthreads();

      real block_result = 0.0;

      for(size_t q = 0; q < BLOCK_SIZE; q++)
      {
         block_result += quads_shared[q].weight * (
            calcLaplaceIntegralGPU(
            quads_shared[q].quad.x, quads_shared[q].quad.y, quads_shared[q].quad.z,
            points_shared[threadIdx.x].x, points_shared[threadIdx.x].y, points_shared[threadIdx.x].z,
            quads_shared[q].normal.x, quads_shared[q].normal.y, quads_shared[q].normal.z));
      }

      point_result += block_result;
      __syncthreads();
   }

   results[point_idx_g] = point_result / fourPI;
}