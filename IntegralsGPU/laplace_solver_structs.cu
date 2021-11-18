#include "laplace_solver_structs.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "laplace_data.cuh"
#include "cuda_helper.h"
#include "laplace_solver_kernels.cuh"
//#include <omp.h>

LaplaceSolverStructs::LaplaceSolverStructs() {};

using namespace laplace_data;

void LaplaceSolverStructs::PrepareData(vector<Vector3>& points, Mesh& mesh, BasisQuadratures& basisQuads)
{
   quadraturesCount = basisQuads.order * mesh.TrianglesCount();
   trianglesCount = mesh.TrianglesCount();
   pointsCount = points.size();
   quadraturesOrder = basisQuads.order;

   // Preparing quadPoints, normals and weights
   quadPoints.resize(quadraturesCount);

   for(size_t t = 0; t < trianglesCount; t++)
   {
      Triangle tr = mesh.GetTriangle(t);

      for(size_t q = 0; q < basisQuads.order; q++)
      {
         int ind = t * basisQuads.order + q;

         quadPoints[ind].quad = tr.PointFromST(basisQuads.x[q], basisQuads.y[q]);
         quadPoints[ind].normal = tr.Normal();
         quadPoints[ind].weight = basisQuads.w[q] * tr.Area();
      }
   }

   // Preparing points
   this->points = vector<Vector3>(points);

   // Preparing results
   results = vector<float>(pointsCount, 0);
}

void LaplaceSolverStructs::CopyToDevice()
{
   // Copying quadPoints
   dev_quadPoints = DevPtr<QuadPoint>(quadPoints.data(), quadraturesCount);

   // Copying points
   dev_points = DevPtr<Vector3>(points.data(), pointsCount);

   // Copying results
   dev_results = DevPtr<float>(pointsCount);
}

vector<float>& LaplaceSolverStructs::SolveCPU()
{
//   size_t p;
//   int n_threads;
//#pragma omp parallel default(none) private(p) shared(trianglesCount, quadraturesOrder, quadPoints, points, results, n_threads)
   {
//#pragma omp master
//      n_threads = omp_get_num_threads();

//#pragma omp for
      for(size_t p = 0; p < pointsCount; p++)
      {
         float point_sum = 0;

         for(size_t t = 0; t < trianglesCount; t++)
         {
            float tringle_sum_1 = 0;
            float tringle_sum_2 = 0;

            for(size_t q = 0; q < quadraturesOrder; q++)
            {
               int idx = t * quadraturesOrder + q;

               tringle_sum_1 += quadPoints[idx].weight * laplaceIntegral1(quadPoints[idx].quad,
                                                                          points[p],
                                                                          quadPoints[idx].normal);

               tringle_sum_2 += quadPoints[idx].weight * laplaceIntegral2(quadPoints[idx].quad,
                                                                          points[p],
                                                                          quadPoints[idx].normal);
            }

            point_sum += tringle_sum_1 - tringle_sum_2;
         }

         results[p] = point_sum / (4.0 * PI);
      }
   }

   return results;
}

void LaplaceSolverStructs::SolveGPU()
{
   switch(algorythmGPU)
   {
      case AlgorythmGPU::Reduction:
      { 
         dim3 dimBlock(THREADS_PER_BLOCK, POINTS_PER_BLOCK);
         dim3 dimGrid(1, pointsCount / POINTS_PER_BLOCK);

         laplace_solver_kernels::SolverKernelStructsReduction<<<
            dimGrid,
            dimBlock,
            THREADS_PER_BLOCK * POINTS_PER_BLOCK * sizeof(float)>>>(
               dev_quadPoints.Get(),
               dev_points.Get(),
               quadraturesCount, dev_results.Get());

         break;
      }
      case AlgorythmGPU::Blocks:
      {
         /*dim3 dimBlock(QUADS_PER_BLOCK, POINTS_PER_BLOCK);
         dim3 dimGrid(1, pointsCount / POINTS_PER_BLOCK);*/

         dim3 dimBlock(POINTS_PER_BLOCK);
         dim3 dimGrid(pointsCount / POINTS_PER_BLOCK);

         laplace_solver_kernels::SolverKernelStructsBlocks<<<
            dimGrid,
            dimBlock >> > (
               dev_quadPoints.Get(), dev_points.Get(),
               trianglesCount, pointsCount, quadraturesOrder, dev_results.Get());

         break;
      }
   }
   
   tryKernelLaunch();
   tryKernelSynchronize();
}

vector<float>& LaplaceSolverStructs::GetResultGPU()
{
   dev_results.CopyToHost(results.data());
   return results;
}