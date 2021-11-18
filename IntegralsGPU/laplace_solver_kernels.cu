#include "laplace_solver_kernels.cuh";
#include "laplace_data.cuh"

typedef unsigned int uint;

//#if (!defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600)
//#else
//__device__ float atomicAdd(float* a, float b) { return b; }
//#endif

using namespace laplace_data;

__global__ void laplace_solver_kernels::SolverKernelArraysReduction(
   const float* quadratures_X,
   const float* quadratures_Y,
   const float* quadratures_Z,
   const float* normals_X,
   const float* normals_Y,
   const float* normals_Z,
   const float* points_X,
   const float* points_Y,
   const float* points_Z,
   const float* weights,
   const float* areas,
   const int trianglesCount,
   const int pointsCount,
   const int quadraturesOrder,
   float* results)
{
   extern __shared__ float shared_array[];
   float* integral = (__shared__ float*)shared_array;

   uint p = blockIdx.x;
   uint t = threadIdx.x;

   if(t < blockDim.x)
      integral[t] = 0;

   while(t < trianglesCount)
   {
      float tringle_sum_1 = 0;
      float tringle_sum_2 = 0;

      for(size_t o = 0; o < quadraturesOrder; o++)
      {
         int ind = t * quadraturesOrder + o;
         tringle_sum_1 += weights[o] * laplaceIntegral1(
            quadratures_X[ind], quadratures_Y[ind], quadratures_Z[ind],
            points_X[p], points_Y[p], points_Z[p],
            normals_X[t], normals_Y[t], normals_Z[t]);

         tringle_sum_2 += weights[o] * laplaceIntegral2(
            quadratures_X[ind], quadratures_Y[ind], quadratures_Z[ind],
            points_X[p], points_Y[p], points_Z[p],
            normals_X[t], normals_Y[t], normals_Z[t]);
      }

      atomicAdd(&integral[t % blockDim.x], (tringle_sum_1 - tringle_sum_2) * areas[t]);
      t += blockDim.x;
   }

   __syncthreads();
   t = threadIdx.x;

   for(unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
   {
      if(t < s)
      {
         integral[t] += integral[t + s];
      }

      __syncthreads();
   }

   __syncthreads();
   t = threadIdx.x;

   if(t == 0)
   {
      results[p] = integral[0] / (4.0 * PI);
   }
}

__global__ void laplace_solver_kernels::SolverKernelVector3sReduction(
   const Vector3* quadPoints,
   const Vector3* normals,
   const Vector3* points,
   const float* weights,
   const float* areas,
   const int trianglesCount,
   const int pointsCount,
   const int quadraturesOrder,
   float* results)
{
   extern __shared__ float shared_array[];
   float* integral = (__shared__ float*)shared_array;

   uint p = blockIdx.x;
   uint t = threadIdx.x;

   if(t < blockDim.x)
      integral[t] = 0;

   while(t < trianglesCount)
   {
      float tringle_sum_1 = 0;
      float tringle_sum_2 = 0;

      for(size_t o = 0; o < quadraturesOrder; o++)
      {
         int ind = t * quadraturesOrder + o;
         tringle_sum_1 += weights[o] * laplaceIntegral1(
            quadPoints[ind], points[p], normals[t]);

         tringle_sum_2 += weights[o] * laplaceIntegral2(
            quadPoints[ind], points[p], normals[t]);
      }

      atomicAdd(&integral[t % blockDim.x], (tringle_sum_1 - tringle_sum_2) * areas[t]);
      t += blockDim.x;
   }

   __syncthreads();
   t = threadIdx.x;

   for(unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
   {
      if(t < s)
      {
         integral[t] += integral[t + s];
      }

      __syncthreads();
   }

   __syncthreads();
   t = threadIdx.x;

   if(t == 0)
   {
      results[p] = integral[0] / (4.0 * PI);
   }
}

__global__ void laplace_solver_kernels::SolverKernelVector3sBlocks(
   const Vector3* quadPoints,
   const Vector3* normals,
   const Vector3* points,
   const float* weights,
   const float* areas,
   const int trianglesCount,
   const int pointsCount,
   const int quadraturesOrder,
   float* results)
{
   // Slice of points in global memory stored in shared memory
   __shared__ Vector3 points_shared[POINTS_PER_BLOCK];

   // Index of point in global and shared memory
   uint point_idx_g = threadIdx.x + POINTS_PER_BLOCK * blockIdx.x;
   uint point_idx_l = threadIdx.x;

   // Loading points into shared memory
   points_shared[point_idx_l] = points[point_idx_g];

   // Result for each point(thread)
   __shared__ float results_shared[POINTS_PER_BLOCK];

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
      __shared__ float weights_shared[QUADS_PER_BLOCK];

      // Normals for current block
      __shared__ Vector3 normals_shared[TR_PER_BLOCK];

      // Areas for current block
      __shared__ float areas_shared[TR_PER_BLOCK];

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
         float triangle_sum_1 = 0.0;
         float triangle_sum_2 = 0.0;

         // Iterating through every quadrature in triangle
         for(size_t q = 0; q < quadraturesOrder; q++)
         {
            const int idx = t * quadraturesOrder + q;

            triangle_sum_1 += weights_shared[idx] * laplaceIntegral1(
               quads_shared[idx], points_shared[point_idx_l], normals_shared[t]);

            triangle_sum_2 += weights_shared[idx] * laplaceIntegral2(
               quads_shared[idx], points_shared[point_idx_l], normals_shared[t]);
         }

         results_shared[point_idx_l] += (triangle_sum_1 - triangle_sum_2) * areas_shared[t];
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
   //__shared__ float results_shared[POINTS_PER_BLOCK];

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
   //   __shared__ float weights_shared[QUADS_PER_BLOCK];

   //   // Normal for current block
   //   __shared__ Vector3 normals_shared;

   //   // Area for current block
   //   __shared__ float areas_shared;

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
   //      float block_sum_1 = 0.0;
   //      float block_sum_2 = 0.0;

   //      // Iterating through every quadrature in block
   //      for(size_t q = 0; q < QUADS_PER_BLOCK; q++)
   //      {
   //         block_sum_1 += weights_shared[q] * laplaceIntegral1(
   //            quads_shared[q], points_shared[point_idx_l], normals_shared);

   //         block_sum_2 += weights_shared[q] * laplaceIntegral2(
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

__global__ void laplace_solver_kernels::SolverKernelStructsReduction(
   const QuadPoint* quadPoints,
   const Vector3* points,
   const int quadraturesCount,
   float* results)
{
   extern __shared__ float shared_array[];

   const uint point_idx_l = threadIdx.y;
   const uint point_idx_g = threadIdx.y + POINTS_PER_BLOCK * blockIdx.y;
   const int offset = THREADS_PER_BLOCK * point_idx_l;

   float* thread_sums = (__shared__ float*)(shared_array + offset);

   const Vector3 point = points[point_idx_g];
   uint quad_idx_l = threadIdx.x;
   thread_sums[quad_idx_l] = 0;

   float block_sum = 0;
   const size_t quadrs_per_thread = (quadraturesCount + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
   const size_t end = (threadIdx.x + 1) * quadrs_per_thread;

   for(size_t quad_idx_g = threadIdx.x * quadrs_per_thread;
              quad_idx_g < end && quad_idx_g < quadraturesCount;
              quad_idx_g++)
   {
      QuadPoint qp = quadPoints[quad_idx_g];

      block_sum += qp.weight * (laplaceIntegral1(qp.quad, point, qp.normal) -
                                laplaceIntegral2(qp.quad, point, qp.normal));
   }

   thread_sums[quad_idx_l] += block_sum;

   /*extern __shared__ float shared_array[];
float* thread_sums = (__shared__ float*)shared_array;

const Vector3 point = points[blockIdx.x];

uint quad_idx_l = threadIdx.x;

thread_sums[quad_idx_l] = 0;

for(size_t quad_idx_g = threadIdx.x * blockDim.x;
           quad_idx_g < (threadIdx.x + 1) * blockDim.x &&
           quad_idx_g < quadraturesOrder * trianglesCount;
           quad_idx_g++)
{

}

while(quad_idx_l < trianglesCount)
{
   float tringle_sum_1 = 0;
   float tringle_sum_2 = 0;

   for(size_t o = 0; o < quadraturesOrder; o++)
   {
       QuadPoint qp = quadPoints[quad_idx_l * quadraturesOrder + o];

      tringle_sum_1 += qp.weight * laplaceIntegral1(
         qp.quad, point, qp.normal);

      tringle_sum_2 += qp.weight * laplaceIntegral2(
         qp.quad, point, qp.normal);
   }

   thread_sums[quad_idx_l % blockDim.x] +=  tringle_sum_1 - tringle_sum_2;
   quad_idx_l += blockDim.x;
}*/

   __syncthreads();

   for(unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
   {
      if(quad_idx_l < s)
      {
         thread_sums[quad_idx_l] += thread_sums[quad_idx_l + s];
      }

      __syncthreads();
   }

   if(quad_idx_l == 0)
   {
      results[point_idx_g] = thread_sums[0] / (4.0 * PI);
   }
}


__global__ void laplace_solver_kernels::SolverKernelStructsBlocks(
   const QuadPoint* quadPoints,
   const Vector3* points,
   const int trianglesCount,
   const int pointsCount,
   const int quadraturesOrder,
   float* results)
{
   // Index of point in global and shared memory
   const uint point_idx_g = threadIdx.x + POINTS_PER_BLOCK * blockIdx.x;
   const uint point_idx_l = threadIdx.x;

   // Loading points into shared memory
   Vector3 point = points[point_idx_g];

   __syncthreads();

   float res = 0.0;

   // Iterating through every block of quadPoints
   for(size_t b = 0; b < (trianglesCount * quadraturesOrder) / QUADS_PER_BLOCK; b++)
   {
      // Slice of quadPoints in global memory stored in shared memory 
      __shared__ QuadPoint quads_shared[QUADS_PER_BLOCK];

      // Index of quadrature to load from global to shared memory
      uint quad_idx_g = threadIdx.x + QUADS_PER_BLOCK * b;
      uint quad_idx_l = threadIdx.x;

      // Loading quadPoints into shared memory
      while(quad_idx_l < QUADS_PER_BLOCK)
      {
         quads_shared[quad_idx_l] = quadPoints[quad_idx_g];
         quad_idx_l += QUADS_PER_BLOCK;
         quad_idx_g += QUADS_PER_BLOCK;
      }

      __syncthreads();

      // Iterating through every quadrature in block
      for(size_t q = 0; q < QUADS_PER_BLOCK; q++)
      {
         QuadPoint qp = quads_shared[q];

         res += qp.weight * laplaceIntegral1(
            qp.quad, point, qp.normal);

         res -= qp.weight * laplaceIntegral2(
            qp.quad, point, qp.normal);
      }

      __syncthreads();
   }

   results[point_idx_g] = res / (4.0 * PI);
}