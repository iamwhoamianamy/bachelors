#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <chrono>

#include "triangle_quadratures.h"
#include "mesh.h"
#include "laplace_solver.h"
#include "cuda_timer.cuh"
#include "laplace_solver_arrays.cuh"

using namespace std;
using namespace triangle_quadratures;
using namespace cuda_utilities;
namespace lscpu = laplace_solver;

int main()
{
   Mesh mesh;

   try
   {
      mesh.InitFromOBJ("../meshes/icosphere_highres.obj");
   }
   catch(Exeption fileExeption)
   {
      cout << fileExeption;
      exit(1);
   }

   QuadPoints qp;

   try
   {
      qp.InitFromTXT("../quadratures/gauss15_xy.txt", "../quadratures/gauss15_w.txt");
   }
   catch(Exeption fileExeption)
   {
      cout << fileExeption;
      exit(1);
   }

   cudaDeviceProp dev_prop;
   int device_count;
   cudaGetDeviceCount(&device_count);

   if(!device_count)
   {
      cout << "No cuda compatable devices found!" << endl;
      exit(2);
   }

   ofstream fout("data.txt");

   for(size_t points_iteration = 0; points_iteration < 10; points_iteration++)
   {
      const int points_count = pow(2, points_iteration + 1);
      vector<double> res;
      vector<Vector3> points(points_count);

      for(size_t i = 0; i < points_count; i++)
      {
         points[i] = { 0.8 / points_count * (i + 1), 0.20, 0.00 };
      }

      LaplaceSolverArrays laplaceSolverArrays;

      // Preparing quadratures and normals
      auto start = std::chrono::steady_clock::now();
      laplaceSolverArrays.PrepareData(points, mesh, qp);
      auto stop = std::chrono::steady_clock::now();
      auto preparation_time = chrono::duration_cast<chrono::microseconds>(stop - start).count() * 1e-6;

      // Solving on CPU
      start = std::chrono::steady_clock::now();
      res = laplaceSolverArrays.SolveCPU();
      stop = std::chrono::steady_clock::now();
      auto cpu_solving_time = chrono::duration_cast<chrono::microseconds>(stop - start).count() * 1e-6;

      // Copying to device
      start = std::chrono::steady_clock::now();
      laplaceSolverArrays.CopyToDevice();
      stop = std::chrono::steady_clock::now();
      auto copying_time = chrono::duration_cast<chrono::microseconds>(stop - start).count() * 1e-6;

      // Solving on GPU
      CudaTimer cudaTimer;
      cudaTimer.Start();
      laplaceSolverArrays.SolveGPU();
      auto gpu_solving_time = cudaTimer.Ellapsed();

      // Getting results from GPU
      start = std::chrono::steady_clock::now();
      res = laplaceSolverArrays.GetResultGPU();
      stop = std::chrono::steady_clock::now();
      auto getting_results_time = chrono::duration_cast<chrono::microseconds>(stop - start).count() * 1e-6;

      cout << setw(30) << "Points count:" << setw(20) << points_count << endl << endl;

      cout << scientific;
      cout << setw(30) << "Preparing quadratures:" << setw(20) << preparation_time << endl;
      cout << setw(30) << "Solving on CPU:" << setw(20) << cpu_solving_time << endl;
      cout << endl;

      cout << setw(30) << "Copying to device:" << setw(20) << copying_time << endl;
      cout << setw(30) << "Solving on GPU:" << setw(20) << gpu_solving_time << endl;
      cout << setw(30) << "Getting results from GPU:" << setw(20) << getting_results_time << endl;
      double total_gpu_time = copying_time + gpu_solving_time + getting_results_time;

      cout << setw(30) << "Total for GPU:" << setw(20) << total_gpu_time << endl;
      cout << endl;

      cout << setw(30) << "GPU speedup:" << setw(20) << cpu_solving_time / total_gpu_time << endl;
      cout << setw(30) << "Total time:" << setw(20) << cpu_solving_time + total_gpu_time + preparation_time << endl;
      cout << endl;

      fout << points_count << "\t" << cpu_solving_time << "\t" << total_gpu_time << endl;
   }

   //// GPU
   //cout << "GPU computation:" << endl;
   //lsgpu::calcIntegralOverMesh(mesh, qp, points, res);

   //for(size_t i = 0; i < points_count; i++)
   //{
   //   cout << "Point: " << scientific << points[i].x << " " << points[i].y << " " << points[i].z << endl;

   //   double true_value = laplace_solver::u(points[i]);
   //   double calc_value = res[i];
   //   double error = abs((true_value - calc_value) / true_value);

   //   cout << "Integral:" << endl;
   //   cout << fixed;
   //   cout << "True value =       " << setw(16) << true_value << endl;
   //   cout << "Calculated value = " << setw(16) << calc_value << endl;
   //   cout << scientific;
   //   cout << "Error            = " << setw(16) << error << endl;
   //}

   //cout << endl << "-----------------------------------------------" << endl << endl;

   // CPU
   /*auto start = std::chrono::steady_clock::now();
   cout << "CPU computation:" << endl;

   lscpu::calcIntegralOverMesh(mesh, qp, points, res);

   auto stop = std::chrono::steady_clock::now();
   auto ellapsed_time_cpu = chrono::duration_cast<chrono::microseconds>(stop - start).count() * 1e-6;
   cout << "Calculation time: " << ellapsed_time_cpu << endl << endl;

   for(size_t i = 0; i < points_count; i++)
   {
      cout << "Point: " << scientific << points[i].x << " " << points[i].y << " " << points[i].z << endl;

      double true_value = laplace_solver::u(points[i]);
      double calc_value = res[i];
      double error = abs((true_value - calc_value) / true_value);

      cout << "Integral:" << endl;
      cout << fixed;
      cout << "True value =       " << setw(16) << true_value << endl;
      cout << "Calculated value = " << setw(16) << calc_value << endl;
      cout << scientific;
      cout << "Error            = " << setw(16) << error << endl;
   }*/

   /*for(size_t i = 0; i < points_count; i++)
   {
      cout << "Point: " << scientific << points[i].x << " " << points[i].y << " " << points[i].z << endl;

      double true_value = laplace_solver::u(points[i]);
      double calc_value = res[i];
      double error = abs((true_value - calc_value) / true_value);

      cout << "Integral:" << endl;
      cout << fixed;
      cout << "True value =       " << setw(16) << true_value << endl;
      cout << "Calculated value = " << setw(16) << calc_value << endl;
      cout << scientific;
      cout << "Error            = " << setw(16) << error << endl;
   }*/

   return 0;
}