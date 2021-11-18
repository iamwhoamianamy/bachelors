#include "laplace_solver_vector3s.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "laplace_data.cuh"
#include "cuda_helper.h"
#include "laplace_solver_kernels.cuh"

//LaplaceSolverVector3s::LaplaceSolverVector3s(AlgorythmGPU alg) : algorythmGPU(alg) {};
LaplaceSolverVector3s::LaplaceSolverVector3s() {};

using namespace laplace_data;

void LaplaceSolverVector3s::PrepareData(vector<Vector3>& points, Mesh& mesh, BasisQuadratures& basisQuads) 
{
   quadraturesCount = basisQuads.order * mesh.TrianglesCount();
   trianglesCount = mesh.TrianglesCount();
   pointsCount = points.size();
   quadraturesOrder = basisQuads.order;

   // Preparing quadPoints
   quadPoints.resize(quadraturesCount);

   for(size_t t = 0; t < trianglesCount; t++)
   {
      Triangle tr = mesh.GetTriangle(t);

      for(size_t o = 0; o < basisQuads.order; o++)
      {
         int ind = t * basisQuads.order + o;
         quadPoints[ind] = tr.PointFromST(basisQuads.x[o], basisQuads.y[o]);
      }
   }

   // Preparing normals
   normals.resize(trianglesCount);

   for(size_t t = 0; t < trianglesCount; t++)
   {
      Triangle tr = mesh.GetTriangle(t);
      normals[t] = tr.Normal();
   }

   // Preparing points
   this->points = vector<Vector3>(points);

   // Preparing weights
   weights = vector<float>(basisQuads.w);

   // Preparing areas
   areas.resize(trianglesCount);

   for(size_t t = 0; t < trianglesCount; t++)
   {
      areas[t] = mesh.GetTriangle(t).Area();
   }

   // Preparing results
   results = vector<float>(pointsCount, 0);
}

void LaplaceSolverVector3s::CopyToDevice() 
{
   // Copying quadPoints
   dev_quadPoints = DevPtr<Vector3>(quadPoints.data(), quadraturesCount);

   // Copying normals
   dev_normals = DevPtr<Vector3>(normals.data(), trianglesCount);

   // Copying points
   dev_points = DevPtr<Vector3>(points.data(), pointsCount);

   // Copying weights
   dev_weights = DevPtr<float>(weights.data(), weights.size());

   // Copying areas
   dev_areas = DevPtr<float>(areas.data(), trianglesCount);

   // Copying results
   dev_results = DevPtr<float>(pointsCount);
}

vector<float>& LaplaceSolverVector3s::SolveCPU()
{
   for(size_t p = 0; p < pointsCount; p++)
   {
      float integral = 0;

      for(size_t t = 0; t < trianglesCount; t++)
      {
         float tringle_sum_1 = 0;
         float tringle_sum_2 = 0;

         for(size_t o = 0; o < quadraturesOrder; o++)
         {
            int ind = t * quadraturesOrder + o;
            tringle_sum_1 += weights[o] * laplaceIntegral1(quadPoints[ind],
                                                           points[p],
                                                           normals[t]);

            tringle_sum_2 += weights[o] * laplaceIntegral2(quadPoints[ind],
                                                           points[p],
                                                           normals[t]);
         }

         integral += (tringle_sum_1 - tringle_sum_2) * areas[t];
      }

      results[p] = integral / (4.0 * PI);
   }

   return results;
}

void LaplaceSolverVector3s::SolveGPU()
{
   switch(algorythmGPU)
   {
      case AlgorythmGPU::Reduction:
      {
         laplace_solver_kernels::SolverKernelVector3sReduction<<<
            pointsCount,
            THREADS_PER_BLOCK,
            THREADS_PER_BLOCK * sizeof(float)>>>(
               dev_quadPoints.Get(),
               dev_normals.Get(),
               dev_points.Get(),
               dev_weights.Get(), dev_areas.Get(),
               trianglesCount, pointsCount, quadraturesOrder, dev_results.Get());

         break;
      }
      case AlgorythmGPU::Blocks:
      {
         /*dim3 dimBlock(QUADS_PER_BLOCK, POINTS_PER_BLOCK);
         dim3 dimGrid(1, pointsCount / POINTS_PER_BLOCK);*/

         dim3 dimBlock(POINTS_PER_BLOCK);
         dim3 dimGrid(pointsCount / POINTS_PER_BLOCK);

         laplace_solver_kernels::SolverKernelVector3sBlocks<<<
            dimGrid,
            dimBlock >>>(
               dev_quadPoints.Get(),
               dev_normals.Get(),
               dev_points.Get(),
               dev_weights.Get(), dev_areas.Get(),
               trianglesCount, pointsCount, quadraturesOrder, dev_results.Get());

         break;
      }
   }

   tryKernelLaunch();
   tryKernelSynchronize();
}

vector<float>& LaplaceSolverVector3s::GetResultGPU()
{
   dev_results.CopyToHost(results.data());
   return results;
}