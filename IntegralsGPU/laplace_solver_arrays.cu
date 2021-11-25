#include "laplace_solver_arrays.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "laplace_data.cuh"
#include "cuda_helper.h"
#include "laplace_solver_kernels.cuh"

LaplaceSolverArrays::LaplaceSolverArrays() {};

using namespace laplace_data;

void LaplaceSolverArrays::PrepareData(const vector<Vector3>& points,
                                      const Mesh& mesh,
                                      const BasisQuadratures& basisQuads)
{
   quadraturesCount = basisQuads.order * mesh.TrianglesCount();
   trianglesCount = mesh.TrianglesCount();
   pointsCount = points.size();
   quadraturesOrder = basisQuads.order;

   // Preparing quadPoints
   quadratures_X.resize(quadraturesCount);
   quadratures_Y.resize(quadraturesCount);
   quadratures_Z.resize(quadraturesCount);

   for(size_t t = 0; t < trianglesCount; t++)
   {
      Triangle tr = mesh.GetTriangle(t);

      for(size_t o = 0; o < basisQuads.order; o++)
      {
         int ind = t * basisQuads.order + o;
         Vector3 v = tr.PointFromST(basisQuads.x[o], basisQuads.y[o]);

         quadratures_X[ind] = v.x;
         quadratures_Y[ind] = v.y;
         quadratures_Z[ind] = v.z;
      }
   }

   // Preparing normals
   normals_X.resize(trianglesCount);
   normals_Y.resize(trianglesCount);
   normals_Z.resize(trianglesCount);

   for(size_t t = 0; t < trianglesCount; t++)
   {
      Triangle tr = mesh.GetTriangle(t);
      Vector3 normal = tr.Normal();

      normals_X[t] = normal.x;
      normals_Y[t] = normal.y;
      normals_Z[t] = normal.z;
   }

   // Preparing points
   points_X.resize(pointsCount);
   points_Y.resize(pointsCount);
   points_Z.resize(pointsCount);

   for(size_t p = 0; p < pointsCount; p++)
   {
      points_X[p] = points[p].x;
      points_Y[p] = points[p].y;
      points_Z[p] = points[p].z;
   }

   // Preparing weights
   weights = vector<real>(basisQuads.w);

   // Preparing areas
   areas.resize(trianglesCount);

   for(size_t t = 0; t < trianglesCount; t++)
   {
      areas[t] = mesh.GetTriangle(t).Area();
   }

   // Preparing results
   results = vector<real>(pointsCount, 0);
}

void LaplaceSolverArrays::CopyToDevice()
{
   // Copying quadPoints
   dev_quadratures_X = DevPtr<real>(quadratures_X.data(), quadraturesCount);
   dev_quadratures_Y = DevPtr<real>(quadratures_Y.data(), quadraturesCount);
   dev_quadratures_Z = DevPtr<real>(quadratures_Z.data(), quadraturesCount);

   // Copying normals
   dev_normals_X = DevPtr<real>(normals_X.data(), trianglesCount);
   dev_normals_Y = DevPtr<real>(normals_Y.data(), trianglesCount);
   dev_normals_Z = DevPtr<real>(normals_Z.data(), trianglesCount);

   // Copying points
   dev_points_X = DevPtr<real>(points_X.data(), pointsCount);
   dev_points_Y = DevPtr<real>(points_Y.data(), pointsCount);
   dev_points_Z = DevPtr<real>(points_Z.data(), pointsCount);

   // Copying weights
   dev_weights = DevPtr<real>(weights.data(), weights.size());

   // Copying areas
   dev_areas = DevPtr<real>(areas.data(), trianglesCount);

   // Copying results
   dev_results = DevPtr<real>(pointsCount);
}

vector<real>& LaplaceSolverArrays::SolveCPU()
{
   for(size_t p = 0; p < pointsCount; p++)
   {
      real integral = 0;

      for(size_t t = 0; t < trianglesCount; t++)
      {
         real tringle_sum_1 = 0;
         real tringle_sum_2 = 0;

         for(size_t o = 0; o < quadraturesOrder; o++)
         {
            int ind = t * quadraturesOrder + o;
            tringle_sum_1 += weights[o] * laplaceIntegral1CPU(quadratures_X[ind], quadratures_Y[ind], quadratures_Z[ind],
                                                              points_X[p], points_Y[p], points_Z[p],
                                                              normals_X[t], normals_Y[t], normals_Z[t]);

            tringle_sum_2 += weights[o] * laplaceIntegral2CPU(quadratures_X[ind], quadratures_Y[ind], quadratures_Z[ind],
                                                              points_X[p], points_Y[p], points_Z[p],
                                                              normals_X[t], normals_Y[t], normals_Z[t]);
         }

         integral += (tringle_sum_1 - tringle_sum_2) * areas[t];
      }

      results[p] = integral / (4.0 * PI);
   }

   return results;
}


void LaplaceSolverArrays::SolveGPU()
{
   laplace_solver_kernels::solverKernelArraysReduction<<<pointsCount, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(real)>>>(
      dev_quadratures_X.Get(), dev_quadratures_Y.Get(), dev_quadratures_Z.Get(),
      dev_normals_X.Get(), dev_normals_Y.Get(), dev_normals_Z.Get(),
      dev_points_X.Get(), dev_points_Y.Get(), dev_points_Z.Get(),
      dev_weights.Get(), dev_areas.Get(),
      trianglesCount, pointsCount, quadraturesOrder, dev_results.Get());

   tryKernelLaunch();
   tryKernelSynchronize();
}

vector<real>& LaplaceSolverArrays::GetResultGPU()
{
   dev_results.CopyToHost(results.data());
   return results;
}

LaplaceSolverArrays::~LaplaceSolverArrays()
{
}
