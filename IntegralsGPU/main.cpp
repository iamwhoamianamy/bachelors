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

using namespace std;
using namespace triangle_quadratures;
using namespace cuda_utilities;
namespace lscpu = laplace_solver;

void printResults(int pointsCount, vector<Vector3>& points, vector<double>& results)
{
   cout << setw(8) << "i" << setw(16) << "pointx" << setw(16) << "pointy" << setw(16) << "pointz";
   cout << setw(16) << "true" << setw(16) << "calc" << setw(16) << "error" << endl;

   for(size_t i = 0; i < pointsCount; i++)
   {
      cout << fixed << setw(8) << i;
      cout << scientific << setw(16) << points[i].x << setw(16) << points[i].y << setw(16) << points[i].z;

      double true_value = laplace_solver::u(points[i]);
      double calc_value = results[i];
      double error = abs((true_value - calc_value) / true_value);

      cout << fixed << setw(16) << true_value << setw(16) << calc_value;
      cout << scientific << setw(16) << error << endl;
   }
}

enum class LaplaceSolvers
{
   Arrays = 0,
   Vector3s = 1
};

void runLaplaceSolverTests(ofstream& fout, Mesh& mesh, QuadPoints& QuadPoints, LaplaceSolvers choose)
{
   for(size_t points_iteration = 12; points_iteration < 13; points_iteration++)
   {
      const int points_count = pow(2, points_iteration);
      //const int points_count = 10;
      vector<double> cpu_results;
      vector<double> gpu_results;
      vector<Vector3> points(points_count);

      for(size_t i = 0; i < points_count; i++)
      {
         points[i] = { 0.4 / points_count * (i + 1), 0.20, 0.00 };
      }

      LaplaceSolver *laplaceSolver;

      if(choose == LaplaceSolvers::Arrays)
         laplaceSolver = new LaplaceSolverArrays();
      else
         laplaceSolver = new LaplaceSolverVector3s();

      // Preparing quadratures and normals
      auto start = std::chrono::steady_clock::now();
      laplaceSolver->PrepareData(points, mesh, QuadPoints);
      auto stop = std::chrono::steady_clock::now();
      auto preparation_time = chrono::duration_cast<chrono::microseconds>(stop - start).count() * 1e-6;

      // Solving on CPU
      start = std::chrono::steady_clock::now();
      //cpu_results = laplaceSolver->SolveCPU();
      stop = std::chrono::steady_clock::now();
      auto cpu_solving_time = chrono::duration_cast<chrono::microseconds>(stop - start).count() * 1e-6;

      // Copying to device
      start = std::chrono::steady_clock::now();
      laplaceSolver->CopyToDevice();
      stop = std::chrono::steady_clock::now();
      auto copying_time = chrono::duration_cast<chrono::microseconds>(stop - start).count() * 1e-6;

      // Solving on GPU
      CudaTimer cudaTimer;
      cudaTimer.Start();
      laplaceSolver->SolveGPU();
      auto gpu_solving_time = cudaTimer.Ellapsed();

      // Getting results from GPU
      start = std::chrono::steady_clock::now();
      gpu_results = laplaceSolver->GetResultGPU();
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

      double speedup_factor = cpu_solving_time / total_gpu_time;
      cout << setw(30) << "GPU speedup:" << setw(20) << speedup_factor << endl;
      cout << setw(30) << "Total time:" << setw(20) << cpu_solving_time + total_gpu_time + preparation_time << endl;
      cout << endl;

      /*cout << endl << "----------------CPU results:---------------" << endl << endl;
      printResults(points_count, points, cpu_results);*/

      cout << endl << "----------------GPU results:---------------" << endl << endl;
      printResults(points_count, points, gpu_results);

      fout << points_count << "\t";
      fout << cpu_solving_time << "\t";
      fout << total_gpu_time << "\t";
      fout << speedup_factor << endl;
   }
}

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

   QuadPoints quad_points;

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
   fout.open("results/laplace_solver_arrays_test_results.txt");
   cout << "Laplace solver with arrays" << endl << endl;
   runLaplaceSolverTests(fout, mesh, quad_points, LaplaceSolvers::Arrays);
   fout.close();

   /*fout.open("results/laplace_solver_vector3s_test_results.txt");
   cout << endl << "Laplace solver with vector3s" << endl << endl;
   runLaplaceSolverTests(fout, mesh, quad_points, LaplaceSolvers::Vector3s);*/

   return 0;
}