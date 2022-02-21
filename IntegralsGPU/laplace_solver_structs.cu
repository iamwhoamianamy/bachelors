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
}

void LaplaceSolverStructs::CopyToDevice()
{
   // Copying quadPoints
   dev_quadPoints = DevPtr<QuadPoint>(quadPoints.data(), quadraturesCount, BLOCK_SIZE);

   // Copying points
   dev_points = DevPtr<Vector3>(points.data(), pointsCount, BLOCK_SIZE);

   // Copying results
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
            calcLaplaceIntegralCPU(quadPoints[q].quad.x, quadPoints[q].quad.y, quadPoints[q].quad.z,
                                   points[p].x, points[p].y, points[p].z,
                                   quadPoints[q].normal.x, quadPoints[q].normal.y, quadPoints[q].normal.z);
      }

      results[p] = point_sum / (fourPI);
   }

   return results;
}

void LaplaceSolverStructs::SolveGPU()
{
   dim3 dimBlock(BLOCK_SIZE);
   dim3 dimGrid(PointsCountPadded() / BLOCK_SIZE);

   laplace_solver_kernels::solverKernelStructs<<<
      dimGrid,
      dimBlock>>>(
         dev_quadPoints.Get(), dev_points.Get(),
         quadraturesCount, dev_results.Get());
   
   tryKernelLaunch();
   tryKernelSynchronize();
}

vector<real>& LaplaceSolverStructs::GetResultGPU()
{
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
