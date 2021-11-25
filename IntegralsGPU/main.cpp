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
#include "cuda_timer.h"
#include "laplace_solver_arrays.h"
#include "laplace_solver_vector3s.h"
#include "laplace_solver_structs.h"
#include "flops_tests.cuh"

using namespace std;
using namespace triangle_quadratures;
using namespace cuda_utilities;

void printResults(int pointsCount, vector<Vector3>& points, vector<real>& results)
{
   cout << setw(8) << "i" << setw(16) << "pointx" << setw(16) << "pointy" << setw(16) << "pointz";
   cout << setw(16) << "true" << setw(16) << "calc" << setw(16) << "error" << endl;

   for(size_t i = 0; i < pointsCount; i++)
   {
      cout << fixed << setw(8) << i;
      cout << scientific << setw(16) << points[i].x << setw(16) << points[i].y << setw(16) << points[i].z;

      real true_value = laplace_data::u(points[i].x, points[i].y, points[i].z);
      real calc_value = results[i];
      real error = abs(true_value - calc_value) / abs(true_value);

      cout << fixed << setw(16) << true_value << setw(16) << calc_value;
      cout << scientific << setw(16) << error << endl;
   }
}

enum class LaplaceSolvers
{
   Arrays = 0,
   Vector3s = 1,
   Structs = 2
};

void runLaplaceSolverTests(ofstream& fout, Mesh& mesh, BasisQuadratures& basisQuads, LaplaceSolvers choose, AlgorythmGPU alg = AlgorythmGPU::Reduction)
{
   for(size_t points_iteration = 9; points_iteration < 10; points_iteration++)
   {
      const size_t points_count = pow(2, 7);
      vector<real> cpu_results;
      vector<real> gpu_results;
      vector<Vector3> points(points_count);

      for(size_t i = 0; i < points_count; i++)
      {
         points[i] = { 0.4f / points_count * (i + 1), 0.20f, 0.00f };
      }

      LaplaceSolver *laplaceSolver;
      
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
            laplaceSolver->algorythmGPU = alg;

            break;
         }
         case LaplaceSolvers::Structs:
         {
            laplaceSolver = new LaplaceSolverStructs();
            laplaceSolver->algorythmGPU = alg;

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
      //cpu_results = laplaceSolver->SolveCPU();
      stop = ::chrono::steady_clock::now();
      auto cpu_solving_time = chrono::duration_cast<chrono::microseconds>(stop - start).count() * 1e-6;

      // Copying to device
      start = ::chrono::steady_clock::now();
      laplaceSolver->CopyToDevice();
      stop = ::chrono::steady_clock::now();
      auto copying_time = chrono::duration_cast<chrono::microseconds>(stop - start).count() * 1e-6;

      // Solving on GPU
      start = ::chrono::steady_clock::now();
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

      if(0)
      {
         cout << endl << "----------------CPU results:---------------" << endl << endl;
         printResults(points_count, points, cpu_results);
      }
      
      if(1)
      {
         cout << endl << "----------------GPU results:---------------" << endl << endl;
         printResults(points_count, points, gpu_results);
      }

      fout << fixed << setprecision(6);

      fout << setw(16) << points_count << " ";
      fout << setw(16) << cpu_solving_time << " ";
      fout << setw(16) << gpu_solving_time << " ";
      fout << setw(16) << speedup_factor << endl;

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
   /*fout.open("results/laplace_solver_arrays_test_results.txt");
   cout << "Laplace solver with arrays" << endl << endl;
   runLaplaceSolverTests(fout, mesh, quad_points, LaplaceSolvers::Arrays);
   fout.close();*/

   /*fout.open("results/laplace_solver_vector3s_reduction_test_results.txt");
   cout << endl << "Laplace solver with vector3s, alg = reduction" << endl << endl;
   runLaplaceSolverTests(fout, mesh, quad_points, LaplaceSolvers::Vector3s, AlgorythmGPU::Reduction);
   fout.close();*/

   //fout.open("results/laplace_solver_vector3s_blocks_test_results.txt");
   //cout << endl << "Laplace solver with vector3s, alg = blocks" << endl << endl;
   //runLaplaceSolverTests(fout, mesh, quad_points, LaplaceSolvers::Vector3s, AlgorythmGPU::Blocks);
   //fout.close();

   //fout.open("results/laplace_solver_structs_reduction_test_results.txt");
   //cout << endl << "Laplace solver with structs, alg = reduction" << endl << endl;
   //runLaplaceSolverTests(fout, mesh, quad_points, LaplaceSolvers::Structs, AlgorythmGPU::Reduction);
   //fout.close();

   /*fout.open("results/laplace_solver_structs_blocks_test_results.txt");
   cout << endl << "Laplace solver with structs, alg = blocks" << endl << endl;
   runLaplaceSolverTests(fout, mesh, quad_points, LaplaceSolvers::Structs, AlgorythmGPU::Blocks);
   fout.close();*/

   fout.open("results/laplace_solver_structs_grid_test_results.txt");
   cout << endl << "Laplace solver with structs, alg = grid" << endl << endl;
   runLaplaceSolverTests(fout, mesh, quad_points, LaplaceSolvers::Structs, AlgorythmGPU::Grid);
   fout.close();
}

int main()
{
   runLaplaceSolverTests();
   
   //addMatricesTest();

   return 0;
}