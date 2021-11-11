#include "laplace_solver_vector3s.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "laplace_data.cuh"
#include "cuda_helper.h"
#include "laplace_solver_kernels.cuh"

LaplaceSolverVector3s::LaplaceSolverVector3s() {};

using namespace laplace_data;

void LaplaceSolverVector3s::PrepareData(vector<Vector3>& points, Mesh& mesh, QuadPoints& quadPoints) 
{
   quadraturesCount = quadPoints.order * mesh.TrianglesCount();
   trianglesCount = mesh.TrianglesCount();
   pointsCount = points.size();
   quadPointsOrder = quadPoints.order;

   // Preparing quadratures
   quadratures.resize(quadraturesCount);

   for(size_t t = 0; t < trianglesCount; t++)
   {
      Triangle tr = mesh.GetTriangle(t);

      for(size_t o = 0; o < quadPoints.order; o++)
      {
         int ind = t * quadPoints.order + o;
         quadratures[ind] = tr.PointFromST(quadPoints.x[o], quadPoints.y[o]);
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

void LaplaceSolverVector3s::CopyToDevice() 
{
   // Copying quadratures
   dev_quadratures = DevPtr<Vector3>(quadratures.data(), quadraturesCount);

   // Copying normals
   dev_normals = DevPtr<Vector3>(normals.data(), trianglesCount);

   // Copying points
   dev_points = DevPtr<Vector3>(points.data(), pointsCount);

   // Copying weights
   dev_weights = DevPtr<double>(weights.data(), weights.size());

   // Copying areas
   dev_areas = DevPtr<double>(areas.data(), trianglesCount);

   // Copying result
   dev_result = DevPtr<double>(pointsCount);
}

vector<double>& LaplaceSolverVector3s::SolveCPU()
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
            tringle_sum_1 += weights[o] * laplaceIntegral1(quadratures[ind],
                                                           points[p],
                                                           normals[t]);

            tringle_sum_2 += weights[o] * laplaceIntegral2(quadratures[ind],
                                                           points[p],
                                                           normals[t]);
         }

         integral += (tringle_sum_1 - tringle_sum_2) * areas[t];
      }

      result[p] = integral / (4.0 * PI);
   }

   return result;
}

void LaplaceSolverVector3s::SolveGPU()
{
   laplace_solver_kernels::SolverKernelVector3s<<<pointsCount, 512, 512 * sizeof(double)>>>(
      dev_quadratures.Get(),
      dev_normals.Get(),
      dev_points.Get(),
      dev_weights.Get(), dev_areas.Get(),
      trianglesCount, pointsCount, quadPointsOrder, dev_result.Get());

   tryKernelLaunch();
   tryKernelSynchronize();
}

vector<double>& LaplaceSolverVector3s::GetResultGPU()
{
   dev_result.CopyToHost(result.data());
   return result;
}