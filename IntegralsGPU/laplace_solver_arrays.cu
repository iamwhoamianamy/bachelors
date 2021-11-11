#include "laplace_solver_arrays.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "laplace_data.cuh"
#include "cuda_helper.h"
#include "laplace_solver_kernels.cuh"

LaplaceSolverArrays::LaplaceSolverArrays() {};

using namespace laplace_data;

void LaplaceSolverArrays::PrepareData(vector<Vector3>& points, Mesh& mesh, QuadPoints& quadPoints)
{
   quadraturesCount = quadPoints.order * mesh.TrianglesCount();
   trianglesCount = mesh.TrianglesCount();
   pointsCount = points.size();
   quadPointsOrder = quadPoints.order;

   // Preparing quadratures
   quadratures_X.resize(quadraturesCount);
   quadratures_Y.resize(quadraturesCount);
   quadratures_Z.resize(quadraturesCount);

   for(size_t t = 0; t < trianglesCount; t++)
   {
      Triangle tr = mesh.GetTriangle(t);

      for(size_t o = 0; o < quadPoints.order; o++)
      {
         int ind = t * quadPoints.order + o;
         Vector3 v = tr.PointFromST(quadPoints.x[o], quadPoints.y[o]);

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
   weights = vector<double>(quadPoints.w);

   // Preparing areas
   areas.resize(trianglesCount);

   for(size_t t = 0; t < trianglesCount; t++)
   {
      areas[t] = mesh.GetTriangle(t).Area();
   }

   // Preparing result
   result = vector<double>(pointsCount, 0);
}

void LaplaceSolverArrays::CopyToDevice()
{
   // Copying quadratures
   dev_quadratures_X = DevPtr<double>(quadratures_X.data(), quadraturesCount);
   dev_quadratures_Y = DevPtr<double>(quadratures_Y.data(), quadraturesCount);
   dev_quadratures_Z = DevPtr<double>(quadratures_Z.data(), quadraturesCount);

   // Copying normals
   dev_normals_X = DevPtr<double>(normals_X.data(), trianglesCount);
   dev_normals_Y = DevPtr<double>(normals_Y.data(), trianglesCount);
   dev_normals_Z = DevPtr<double>(normals_Z.data(), trianglesCount);

   // Copying points
   dev_points_X = DevPtr<double>(points_X.data(), pointsCount);
   dev_points_Y = DevPtr<double>(points_Y.data(), pointsCount);
   dev_points_Z = DevPtr<double>(points_Z.data(), pointsCount);

   // Copying weights
   dev_weights = DevPtr<double>(weights.data(), weights.size());

   // Copying areas
   dev_areas = DevPtr<double>(areas.data(), trianglesCount);

   // Copying result
   dev_result = DevPtr<double>(pointsCount);
}

vector<double>& LaplaceSolverArrays::SolveCPU()
{
   for(size_t p = 0; p < pointsCount; p++)
   {
      double integral = 0;

      for(size_t t = 0; t < trianglesCount; t++)
      {
         double tringle_sum_1 = 0;
         double tringle_sum_2 = 0;

         for(size_t o = 0; o < quadPointsOrder; o++)
         {
            int ind = t * quadPointsOrder + o;
            tringle_sum_1 += weights[o] * laplaceIntegral1(quadratures_X[ind], quadratures_Y[ind], quadratures_Z[ind],
                                                           points_X[p], points_Y[p], points_Z[p],
                                                           normals_X[t], normals_Y[t], normals_Z[t]);

            tringle_sum_2 += weights[o] * laplaceIntegral2(quadratures_X[ind], quadratures_Y[ind], quadratures_Z[ind],
                                                           points_X[p], points_Y[p], points_Z[p],
                                                           normals_X[t], normals_Y[t], normals_Z[t]);
         }

         integral += (tringle_sum_1 - tringle_sum_2) * areas[t];
      }

      result[p] = integral / (4.0 * PI);
   }

   return result;
}


void LaplaceSolverArrays::SolveGPU()
{
   laplace_solver_kernels::SolverKernelArrays<<<pointsCount, threads_per_block, threads_per_block * sizeof(double)>>>(
      dev_quadratures_X.Get(), dev_quadratures_Y.Get(), dev_quadratures_Z.Get(),
      dev_normals_X.Get(), dev_normals_Y.Get(), dev_normals_Z.Get(),
      dev_points_X.Get(), dev_points_Y.Get(), dev_points_Z.Get(),
      dev_weights.Get(), dev_areas.Get(),
      trianglesCount, pointsCount, quadPointsOrder, dev_result.Get());

   tryKernelLaunch();
   tryKernelSynchronize();
}

vector<double>& LaplaceSolverArrays::GetResultGPU()
{
   dev_result.CopyToHost(result.data());
   return result;
}