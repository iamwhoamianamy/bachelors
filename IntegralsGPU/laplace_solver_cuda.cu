#include <iostream>
#include <chrono>

#include "laplace_solver_cuda.cuh"
#include "dev_ptr.cuh"
#include "cuda_timer.cuh"
#include "cuda_helper.cuh"

//#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
//#else
//__device__ double atomicAdd(double* a, double b) { return b; }
//#endif

using namespace cuda_utilities;

__device__ const double PI = 3.14159265359;

__device__ double laplace_solver_cuda::u(double x, double y, double z)
{
   return 2 * x * x -y * y - z * z;
}

__device__ double laplace_solver_cuda::gradUX(double x, double y, double z)
{
   return 4 * x;
}

__device__ double laplace_solver_cuda::gradUY(double x, double y, double z)
{
   return -2 * y;
}

__device__ double laplace_solver_cuda::gradUZ(double x, double y, double z)
{
   return -2 * z;
}

__device__ double laplace_solver_cuda::lengthBetween(double x1, double y1, double z1,
                                                     double x2, double y2, double z2)
{
   return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
}

__device__ double laplace_solver_cuda::laplaceIntegral1(double qx, double qy, double qz,
                                                          double px, double py, double pz,
                                                          double nx, double ny, double nz)
{
   double gradx = gradUX(qx, qy, qz);
   double grady = gradUY(qx, qy, qz);
   double gradz = gradUZ(qx, qy, qz);

   double dudnx = gradx * nx;
   double dudny = grady * ny;
   double dudnz = gradz * nz;

   return (dudnx + dudny + dudnz) / lengthBetween(qx, qy, qz, px, py, pz);
}

__device__ double laplace_solver_cuda::laplaceIntegral2(double qx, double qy, double qz,
                                                        double px, double py, double pz,
                                                        double nx, double ny, double nz)
{
   double l = lengthBetween(qx, qy, qz, px, py, pz);

   double rx = nx * (px - qx);
   double ry = ny * (py - qy);
   double rz = nz * (pz - qz);

   return (rx + ry + rz) * u(qx, qy, qz) / (l * l * l);
}

void laplace_solver_cuda::calcIntegralOverMesh(const Mesh& mesh,
                                               const QuadPoints& qp,
                                               const vector<Vector3>& points,
                                               vector<double>& result)
{
   auto start = chrono::steady_clock::now();

   const int quadraturesCount = qp.order * mesh.TrianglesCount();

   // Preparing quadratures
   vector<double> quadraturesX(quadraturesCount);
   vector<double> quadraturesY(quadraturesCount);
   vector<double> quadraturesZ(quadraturesCount);

   for(size_t t = 0; t < mesh.TrianglesCount(); t++)
   {
      Triangle tr = mesh.GetTriangle(t);

      for(size_t o = 0; o < qp.order; o++)
      {
         int ind = t * qp.order + o;
         Vector3 v = tr.PointFromST(qp.x[o], qp.y[o]);

         quadraturesX[ind] = v.x;
         quadraturesY[ind] = v.y;
         quadraturesZ[ind] = v.z;
      }
   }

   DevPtr<double> dev_quadraturesX(quadraturesX.data(), quadraturesCount);
   DevPtr<double> dev_quadraturesY(quadraturesY.data(), quadraturesCount);
   DevPtr<double> dev_quadraturesZ(quadraturesZ.data(), quadraturesCount);

   // Preparing normals
   vector<double> normalsX(mesh.TrianglesCount());
   vector<double> normalsY(mesh.TrianglesCount());
   vector<double> normalsZ(mesh.TrianglesCount());

   for(size_t t = 0; t < mesh.TrianglesCount(); t++)
   {
      Triangle tr = mesh.GetTriangle(t);
      Vector3 normal = tr.Normal();
      
      normalsX[t] = normal.x;
      normalsY[t] = normal.y;
      normalsZ[t] = normal.z;
   }

   DevPtr<double> dev_normalsX(normalsX.data(), mesh.TrianglesCount());
   DevPtr<double> dev_normalsY(normalsY.data(), mesh.TrianglesCount());
   DevPtr<double> dev_normalsZ(normalsZ.data(), mesh.TrianglesCount());

   // Preparing points
   vector<double> pointsX(points.size());
   vector<double> pointsY(points.size());
   vector<double> pointsZ(points.size());

   for(size_t p = 0; p < points.size(); p++)
   {
      pointsX[p] = points[p].x;
      pointsY[p] = points[p].y;
      pointsZ[p] = points[p].z;
   }

   DevPtr<double> dev_points_X(pointsX.data(), points.size());
   DevPtr<double> dev_points_Y(pointsY.data(), points.size());
   DevPtr<double> dev_points_Z(pointsZ.data(), points.size());

   // Preparing weights
   DevPtr<double> dev_weights(qp.w.data(), qp.order);

   // Preparing areas
   vector<double> areas(mesh.TrianglesCount());

   for(size_t t = 0; t < mesh.TrianglesCount(); t++)
   {
      areas[t] = mesh.GetTriangle(t).Area();
   }

   DevPtr<double> dev_areas(areas.data(), mesh.TrianglesCount());

   // Preparing result
   result = vector<double>(points.size(), 0);
   DevPtr<double> dev_result(points.size());

   auto stop = std::chrono::steady_clock::now();
   auto ellapsed_time_cpu = chrono::duration_cast<chrono::microseconds>(stop - start).count() * 1e-6;
   cout << "Preparation duration: " << scientific << ellapsed_time_cpu << endl;
   
   //dim3 grid_blocks(points.size(), mesh.TrianglesCount());
   //dim3 grid_threads(512, 512);

   CudaTimer cuda_timer;
   cuda_timer.Start();

   /*calcIntegralOverMeshArrays<<<points.size(), 512, 512 * sizeof(double)>>>(
      dev_quadraturesX.Get(), dev_quadraturesY.Get(), dev_quadraturesZ.Get(),
      dev_normalsX.Get(), dev_normalsY.Get(), dev_normalsZ.Get(),
      dev_points_X.Get(), dev_points_Y.Get(), dev_points_Z.Get(),
      dev_weights.Get(), dev_areas.Get(),
      mesh.TrianglesCount(), points.size(), qp.order, dev_result.Get());*/

   tryKernelLaunch();
   tryKernelSynchronize();

   auto ellapsed_time_gpu = cuda_timer.Ellapsed();
   cout << "Kernel runtime: " << ellapsed_time_gpu << endl;
   cout << "Calculation time: " << ellapsed_time_cpu + ellapsed_time_gpu  << endl;

   dev_result.CopyToHost(result.data());
}

//__global__ void laplace_solver_cuda::calcIntegralOverMeshArrays(const double* quadraturesX,
//                                                                const double* quadraturesY,
//                                                                const double* quadraturesZ,
//                                                                const double* normalsX,
//                                                                const double* normalsY,
//                                                                const double* normalsZ,
//                                                                const double* pointsX,
//                                                                const double* pointsY,
//                                                                const double* pointsZ,
//                                                                const double* weights,
//                                                                const double* areas,
//                                                                const int trianglesCount,
//                                                                const int pointsCount,
//                                                                const int quadratureOrder,
//                                                                double* result)
//{
//   uint p = blockIdx.x;
//   extern __shared__ double shared_array[];
//   double* integral = (__shared__ double*)shared_array;
//   uint t = threadIdx.x;
//
//   while(t < trianglesCount)
//   {
//      double tringle_sum_1 = 0;
//      double tringle_sum_2 = 0;
//
//      for(size_t o = 0; o < quadratureOrder; o++)
//      {
//         int ind = t * quadratureOrder + o;
//         tringle_sum_1 += weights[o] * laplaceIntegral1(quadraturesX[ind], quadraturesY[ind], quadraturesZ[ind],
//                                                         pointsX[p], pointsY[p], pointsZ[p],
//                                                         normalsX[t], normalsY[t], normalsZ[t]);
//
//         tringle_sum_2 += weights[o] * laplaceIntegral2(quadraturesX[ind], quadraturesY[ind], quadraturesZ[ind],
//                                                         pointsX[p], pointsY[p], pointsZ[p],
//                                                         normalsX[t], normalsY[t], normalsZ[t]);
//      }
//
//      atomicAdd(&integral[t % blockDim.x], (tringle_sum_1 - tringle_sum_2) * areas[t]);
//      t += blockDim.x;
//   }
//
//   __syncthreads();
//   t = threadIdx.x;
//
//   for(unsigned int s = 1; s < blockDim.x; s *= 2)
//   {
//      if(t % (2 * s) == 0)
//      {
//         if(t + s < trianglesCount)
//         {
//            integral[t] += integral[t + s];
//         }
//      }
//
//      __syncthreads();
//   }
//
//   __syncthreads();
//   t = threadIdx.x;
//
//   if(t == 0)
//   {
//      result[p] = integral[0] / (4.0 * PI);
//   }
//}
