#include "laplace_solver_structs.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "laplace_data.cuh"
#include "cuda_helper.h"
#include "laplace_solver_kernels.cuh"
//#include <omp.h>

LaplaceSolverStructs::LaplaceSolverStructs() {};

using namespace laplace_data;

void LaplaceSolverStructs::PrepareData(const vector<Vector3>& points,
                                       const Mesh& mesh,
                                       const BasisQuadratures& basisQuads)
{
   // Initializing constants
   quadraturesCount = basisQuads.order * mesh.TrianglesCount();
   pointsCount = points.size();

   // Preparing quadPoints, normals and weights
   quadPoints.resize(quadraturesCount);

   for(size_t t = 0; t < mesh.TrianglesCount(); t++)
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
   results = vector<real>(pointsCount, 0.0f);

   if(algorythmGPU == AlgorythmGPU::Grid)
   {
      matrixWidth = nextDevisible(quadraturesCount, BLOCK_SIZE) / BLOCK_SIZE;
      resultsMatrix = vector<real>(PointsCountPadded() * matrixWidth, 0.0f);
   }
}

void LaplaceSolverStructs::CopyToDevice()
{
   // Copying quadPoints
   dev_quadPoints = DevPtr<QuadPoint>(quadPoints.data(), quadraturesCount, BLOCK_SIZE);

   // Copying points
   dev_points = DevPtr<Vector3>(points.data(), pointsCount, BLOCK_SIZE);

   // Copying results
   if(algorythmGPU == AlgorythmGPU::Grid)
      dev_resultsMatrix = DevPtr<real>(resultsMatrix.size());
   else
      dev_results = DevPtr<real>(pointsCount, BLOCK_SIZE);

}

vector<real>& LaplaceSolverStructs::SolveCPU()
{
   for(size_t p = 0; p < pointsCount; p++)
   {
      real point_sum = 0;

      for(size_t q = 0; q < quadraturesCount; q++)
      {
         point_sum += quadPoints[q].weight *
            (laplaceIntegral1CPU(quadPoints[q].quad, points[p], quadPoints[q].normal) -
             laplaceIntegral2CPU(quadPoints[q].quad, points[p], quadPoints[q].normal));
      }

      results[p] = point_sum / (4.0 * PI);
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

         laplace_solver_kernels::solverKernelStructsReduction<<<
            dimGrid,
            dimBlock,
            THREADS_PER_BLOCK * POINTS_PER_BLOCK * sizeof(real)>>>(
               dev_quadPoints.Get(),
               dev_points.Get(),
               quadraturesCount, dev_results.Get());

         break;
      }
      case AlgorythmGPU::Blocks:
      {
         dim3 dimBlock(BLOCK_SIZE);
         dim3 dimGrid(PointsCountPadded() / BLOCK_SIZE);

         laplace_solver_kernels::solverKernelStructsBlocks<<<
            dimGrid,
            dimBlock >> > (
               dev_quadPoints.Get(), dev_points.Get(),
               quadraturesCount, dev_results.Get());

         break;
      }
      case AlgorythmGPU::Grid:
      {
         dim3 dimBlock(1, BLOCK_SIZE);
         dim3 dimGrid(matrixWidth, PointsCountPadded() / BLOCK_SIZE);

         laplace_solver_kernels::solverKernelStructsGrid<<<
            dimGrid,
            dimBlock>>>(
               dev_quadPoints.Get(), dev_points.Get(),
               matrixWidth, dev_resultsMatrix.Get());

         break;
      }
   }
   
   tryKernelLaunch();
   tryKernelSynchronize();
}

vector<real>& LaplaceSolverStructs::GetResultGPU()
{
   if(algorythmGPU == AlgorythmGPU::Grid)
   {
      dev_resultsMatrix.CopyToHost(resultsMatrix.data());

      for(size_t i = 0; i < pointsCount; i++)
      {
         results[i] = 0.0f;

         for(size_t j = 0; j < matrixWidth; j++)
         {
            results[i] += resultsMatrix[i * matrixWidth + j];
         }

         results[i] /= (4.0 * PI);
      }
   }
   else
      dev_results.CopyToHost(results.data());

   return results;
}

size_t LaplaceSolverStructs::PointsCountPadded() const
{
   return nextDevisible(pointsCount, BLOCK_SIZE);
}

#include <iostream>

LaplaceSolverStructs::~LaplaceSolverStructs()
{
   cout << "LaplaceSolverStruct was destroyed!" << endl;
}
