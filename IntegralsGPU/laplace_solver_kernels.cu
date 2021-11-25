#include "laplace_solver_kernels.cuh";
#include "laplace_data.cuh"

typedef unsigned int uint;

#if (!defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600)
#else
__device__ real atomicAdd(double* a, double b) { return b; }
#endif

using namespace laplace_data;

__global__ void laplace_solver_kernels::solverKernelArraysReduction(
   const real* quadratures_X,
   const real* quadratures_Y,
   const real* quadratures_Z,
   const real* normals_X,
   const real* normals_Y,
   const real* normals_Z,
   const real* points_X,
   const real* points_Y,
   const real* points_Z,
   const real* weights,
   const real* areas,
   const int trianglesCount,
   const int pointsCount,
   const int quadraturesOrder,
   real* results)
{
   extern __shared__ real shared_array[];
   real* integral = (__shared__ real*)shared_array;

   uint p = blockIdx.x;
   uint t = threadIdx.x;

   if(t < blockDim.x)
      integral[t] = 0;

   while(t < trianglesCount)
   {
      real tringle_sum_1 = 0;
      real tringle_sum_2 = 0;

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

__global__ void laplace_solver_kernels::solverKernelVector3sReduction(
   const Vector3* quadPoints,
   const Vector3* normals,
   const Vector3* points,
   const real* weights,
   const real* areas,
   const int trianglesCount,
   const int pointsCount,
   const int quadraturesOrder,
   real* results)
{
   extern __shared__ real shared_array[];
   real* integral = (__shared__ real*)shared_array;

   uint p = blockIdx.x;
   uint t = threadIdx.x;

   if(t < blockDim.x)
      integral[t] = 0;

   while(t < trianglesCount)
   {
      real tringle_sum_1 = 0;
      real tringle_sum_2 = 0;

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

__global__ void laplace_solver_kernels::solverKernelVector3sBlocks(
   const Vector3* quadPoints,
   const Vector3* normals,
   const Vector3* points,
   const real* weights,
   const real* areas,
   const int trianglesCount,
   const int pointsCount,
   const int quadraturesOrder,
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
         real triangle_sum_1 = 0.0;
         real triangle_sum_2 = 0.0;

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

__global__ void laplace_solver_kernels::solverKernelStructsReduction(
   const QuadPoint* quadPoints,
   const Vector3* points,
   const int quadraturesCount,
   real* results)
{
   extern __shared__ real shared_array[];

   const uint point_idx_l = threadIdx.y;
   const uint point_idx_g = threadIdx.y + POINTS_PER_BLOCK * blockIdx.y;
   const int offset = THREADS_PER_BLOCK * point_idx_l;

   real* thread_sums = (__shared__ real*)(shared_array + offset);

   const Vector3 point = points[point_idx_g];
   uint quad_idx_l = threadIdx.x;
   thread_sums[quad_idx_l] = 0;

   real block_sum = 0;
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

__global__ void laplace_solver_kernels::solverKernelStructsBlocks(
   const QuadPoint* quadPoints,
   const Vector3* points,
   const int quadraturesCount,
   real* results)
{
   // Index of point in global and shared memory
   const uint point_idx_g = threadIdx.x + BLOCK_SIZE * blockIdx.x;
   const uint point_idx_l = threadIdx.x;

   // Loading points into shared memory
   Vector3 point = points[point_idx_g];

   __syncthreads();

   real result = 0.0;

   // Iterating through every block of quadPoints
   for(size_t b = 0; b < quadraturesCount / BLOCK_SIZE; b++)
   {
      // Slice of quadPoints in global memory stored in shared memory 
      __shared__ QuadPoint quads_shared[BLOCK_SIZE];

      // Index of quadrature to load from global to shared memory
      const uint quad_idx_g = threadIdx.x + BLOCK_SIZE * b;
      const uint quad_idx_l = threadIdx.x;

      // Loading quadPoints into shared memory
      quads_shared[quad_idx_l] = quadPoints[quad_idx_g];

      __syncthreads();

      // Iterating through every quadrature in block
      for(size_t q = 0; q < BLOCK_SIZE; q++)
      {
         QuadPoint qp = quads_shared[q];

         result += qp.weight * (
            laplaceIntegral1( qp.quad, point, qp.normal) -
            laplaceIntegral2( qp.quad, point, qp.normal));
      }

      __syncthreads();
   }

   results[point_idx_g] = result / (4.0 * PI);
}

__global__ void laplace_solver_kernels::solverKernelStructsGrid(
   const QuadPoint* quadPoints,
   const Vector3* points,
   const int matrixWidth,
   real* resultsMatrix)
{
   __shared__ QuadPoint quads_shared[BLOCK_SIZE];

   const uint point_idx_g = threadIdx.y + blockDim.y * blockIdx.y;

   const Vector3 point = points[point_idx_g];
   quads_shared[threadIdx.y] = quadPoints[threadIdx.y + blockDim.y * blockIdx.x];

   __syncthreads();

   real result = 0.0;

   for(size_t q = 0; q < BLOCK_SIZE; q++)
   {
      result += quads_shared[q].weight * (
            laplaceIntegral1(quads_shared[q].quad, point, quads_shared[q].normal) -
            laplaceIntegral2(quads_shared[q].quad, point, quads_shared[q].normal));
   }

   resultsMatrix[point_idx_g * matrixWidth + blockIdx.x] = result;
}

__global__ void laplace_solver_kernels::AddMatrices(const real* a, const real* b, real* c)
{
   uint i = threadIdx.y + blockDim.y * blockIdx.y;
   uint j = threadIdx.x + blockDim.x * blockIdx.x;

   c[i * MATRIX_WIDTH + j] = a[i * MATRIX_WIDTH + j] + b[i * MATRIX_WIDTH + j];
}

__global__ void laplace_solver_kernels::AddMatricesShared(const real* a, const real* b, real* c)
{
   __shared__ real a_sub[BLOCK_SIZE][BLOCK_SIZE];
   __shared__ real b_sub[BLOCK_SIZE][BLOCK_SIZE];
   __shared__ real c_sub[BLOCK_SIZE][BLOCK_SIZE];

   uint i_g = threadIdx.y + blockDim.y * blockIdx.y;
   uint j_g = threadIdx.x + blockDim.x * blockIdx.x;

   uint i_l = threadIdx.y;
   uint j_l = threadIdx.x;

   a_sub[i_l][j_l] = a[i_g * MATRIX_WIDTH + j_g];
   b_sub[i_l][j_l] = b[i_g * MATRIX_WIDTH + j_g];

   __syncthreads();

   real k = 0.0;

   while(k < 200000.0)
   {
      c_sub[i_l][j_l] = a_sub[i_l][j_l] + c_sub[i_l][j_l];
      k += 1.0;
   }
   
   __syncthreads();

   c[i_g * MATRIX_WIDTH + j_g] = c_sub[i_l][j_l];
}