﻿#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <chrono>

#include "triangle_quadratures.h"
#include "mesh.h"
#include "laplace_data.cuh"
#include "cuda_timer.h"
#include "laplace_solver_arrays.h"
#include "laplace_solver_vector3s.h"
#include "laplace_solver_structs.h"

using namespace std;
using namespace triangle_quadratures;
using namespace cuda_utilities;

real calcAverageError(int pointsCount, vector<Vector3>& points, vector<real>& results)
{
   real averageError = 0;

   for(size_t i = 0; i < pointsCount; i++)
   {
      real true_value = laplace_data::u(points[i].x, points[i].y, points[i].z);
      real calc_value = results[i];
      real error = abs(true_value - calc_value) / abs(true_value);

      averageError += error;
   }

   return averageError / pointsCount;
}

void printResults(int pointsCount, vector<Vector3>& points, vector<real>& results)
{
   cout << setw(6) << "i" << setw(14) << "pointx" << setw(14) << "pointy" << setw(14) << "pointz";
   cout << setw(14) << "true" << setw(14) << "calc" << setw(14) << "error" << endl;

   for(size_t i = 0; i < pointsCount; i++)
   {
      cout << fixed << setw(6) << i;
      cout << scientific << setw(14) << points[i].x << setw(14) << points[i].y << setw(14) << points[i].z;

      real true_value = laplace_data::u(points[i].x, points[i].y, points[i].z);
      real calc_value = results[i];
      real error = abs(true_value - calc_value) / abs(true_value);

      cout << fixed << setw(14) << true_value << setw(14) << calc_value;
      cout << scientific << setw(14) << error << endl;
   }
}

enum class LaplaceSolvers
{
   Arrays = 0,
   Vector3s = 1,
   Structs = 2
};

enum class Devices
{
   CPU = 0,
   GPU = 1
};

void runLaplaceSolverTests(
   ofstream& fout,
   Mesh& mesh,
   BasisQuadratures& basisQuads,
   LaplaceSolvers choose,
   Devices device,
   bool doPrintResults)
{
   for(size_t points_iteration = 10; points_iteration < 11; points_iteration++)
   {
      const size_t points_count = pow(2, points_iteration);
      vector<real> cpu_results;
      vector<real> gpu_results;
      vector<Vector3> points(points_count);

      /*for(size_t i = 0; i < points_count; i++)
      {
         points[i] = { 0.4f / points_count * (i + 1), 0.20f, 0.00f };
      }*/

      for(size_t i = 0; i < points_count; i++)
      {
         points[i] = { 0.9 + 0.1 * (i + 1) / points_count, 0.00f };
      }

      LaplaceSolver* laplaceSolver;

      switch(choose)
      {
         case LaplaceSolvers::Arrays:
         {
            laplaceSolver = new LaplaceSolverArrays();

            break;
         }
         case LaplaceSolvers::Vector3s:
         {
            laplaceSolver = new LaplaceSolverVector3s();

            break;
         }
         case LaplaceSolvers::Structs:
         {
            laplaceSolver = new LaplaceSolverStructs();

            break;
         }
      }

      // Preparing quadPoints and normals
      auto start = ::chrono::steady_clock::now();
      laplaceSolver->PrepareData(points, mesh, basisQuads);
      auto stop = ::chrono::steady_clock::now();
      auto preparation_time = chrono::duration_cast<chrono::microseconds>(stop - start).count() * 1e-6;

      // Solving on CPU
      start = ::chrono::steady_clock::now();
      if(device == Devices::CPU)
         cpu_results = laplaceSolver->SolveCPU();
      stop = ::chrono::steady_clock::now();
      auto cpu_solving_time = chrono::duration_cast<chrono::microseconds>(stop - start).count() * 1e-6;

      // Copying to device
      start = ::chrono::steady_clock::now();
      laplaceSolver->CopyToDevice();
      stop = ::chrono::steady_clock::now();
      auto copying_time = chrono::duration_cast<chrono::microseconds>(stop - start).count() * 1e-6;

      // Solving on GPU
      start = ::chrono::steady_clock::now();
      if(device == Devices::GPU)
         laplaceSolver->SolveGPU();
      stop = ::chrono::steady_clock::now();
      auto gpu_solving_time = chrono::duration_cast<chrono::microseconds>(stop - start).count() * 1e-6;

      // Getting results from GPU
      start = ::chrono::steady_clock::now();
      gpu_results = laplaceSolver->GetResultGPU();
      stop = ::chrono::steady_clock::now();
      auto getting_results_time = chrono::duration_cast<chrono::microseconds>(stop - start).count() * 1e-6;

      cout << "----------------------------------------------------------" << endl;
      cout << setw(30) << "Points count:" << setw(20) << points_count << endl;
      cout << setw(30) << "Triangles count:" << setw(20) << mesh.TrianglesCount() << endl;
      cout << setw(30) << "Quadrature order:" << setw(20) << basisQuads.order << endl << endl;

      cout << scientific;
      cout << setw(30) << "Preparing quadratures:" << setw(20) << preparation_time << endl;
      cout << setw(30) << "Solving on CPU:" << setw(20) << cpu_solving_time << endl;
      cout << endl;

      cout << setw(30) << "Copying to device:" << setw(20) << copying_time << endl;
      cout << setw(30) << "Solving on GPU:" << setw(20) << gpu_solving_time << endl;
      cout << setw(30) << "Getting results from GPU:" << setw(20) << getting_results_time << endl;
      auto total_gpu_time = copying_time + gpu_solving_time + getting_results_time;

      cout << setw(30) << "Total for GPU:" << setw(20) << total_gpu_time << endl;
      cout << endl;

      auto speedup_factor = cpu_solving_time / gpu_solving_time;
      cout << setw(30) << "GPU speedup:" << setw(20) << speedup_factor << endl;
      cout << setw(30) << "Total time:" << setw(20) << cpu_solving_time + total_gpu_time + preparation_time << endl;
      
      if(device == Devices::CPU)
         cout << setw(30) << "Average error CPU:" << setw(20) << calcAverageError(points_count, points, cpu_results);

      if(device == Devices::GPU)
         cout << setw(30) << "Average error GPU:" << setw(20) << calcAverageError(points_count, points, gpu_results);

      if(device == Devices::CPU && doPrintResults)
      {
         cout << endl << "----------------CPU results:---------------" << endl;
         printResults(points_count, points, cpu_results);
      }

      if(device == Devices::GPU && doPrintResults)
      {
         cout << endl << "----------------GPU results:---------------" << endl;
         printResults(points_count, points, gpu_results);
      }

      fout << fixed << setprecision(6);

      //fout << setw(16) << points_count << " ";
      //fout << setw(16) << cpu_solving_time << endl;


      if(device == Devices::CPU)
         fout << cpu_solving_time << endl;

      if(device == Devices::GPU)
         fout << gpu_solving_time << endl;

      //fout << setw(16) << gpu_solving_time << " ";
      //fout << setw(16) << speedup_factor << endl;

      delete laplaceSolver;
   }
}

void runLaplaceSolverTests()
{
   Mesh mesh;

   try
   {
      mesh.InitFromOBJ("../meshes/icospheres/ico20480.obj");
   }
   catch(Exeption fileExeption)
   {
      cout << fileExeption;
      exit(1);
   }

   BasisQuadratures quad_points;

   try
   {
      quad_points.InitFromTXT("../quadratures/gauss15_xy.txt", "../quadratures/gauss15_w.txt");
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

   ofstream fout;

   /*fout.open("results/laplace_solver_structs_results_CPU_double.txt");
   cout << "Laplace solver with structs, alg = blocks" << endl;
   runLaplaceSolverTests(fout, mesh, quad_points, LaplaceSolvers::Structs, Devices::CPU, true);
   fout.close();*/

   //fout.open("results/laplace_solver_structs_results_GPU_float.txt");
   fout.open("results/errors.txt");
   cout << "Laplace solver with structs, alg = blocks" << endl;
   runLaplaceSolverTests(fout, mesh, quad_points, LaplaceSolvers::Structs, Devices::CPU, true);
   fout.close();
}

int main()
{
   runLaplaceSolverTests();
   
   return 0;
}