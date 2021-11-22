#pragma once #include "real.h"
#include <chrono>

#include "flops_tests.cuh"
#include "cuda_timer.h"

void addMatricesTest()
{
   const int MATRIX_LENGTH = MATRIX_WIDTH * MATRIX_HEIGHT;

   double* a = new double[MATRIX_LENGTH];
   double* b = new double[MATRIX_LENGTH];
   double* c = new double[MATRIX_LENGTH];

   int ka = 0;
   int kb = MATRIX_LENGTH;

   for(size_t i = 0; i < MATRIX_HEIGHT; i++)
   {
      for(size_t j = 0; j < MATRIX_WIDTH; j++)
      {
         a[i * MATRIX_WIDTH + j] = ka;
         b[i * MATRIX_WIDTH + j] = kb;

         ka += 0.001f;
         kb -= 0.001f;
      }
   }

   cout << "Done initializing!" << endl;

   DevPtr<double> dev_a(a, MATRIX_LENGTH);
   DevPtr<double> dev_b(b, MATRIX_LENGTH);
   DevPtr<double> dev_c(c, MATRIX_LENGTH);

   cout << "Done copying to GPU!" << endl;

   dim3 gridDim(MATRIX_WIDTH / BLOCK_SIZE, MATRIX_HEIGHT / BLOCK_SIZE);
   dim3 blockDim(BLOCK_SIZE, BLOCK_SIZE);

   auto start = std::chrono::steady_clock::now();
   //AddMatricesShared<<<gridDim, blockDim >>>(dev_a.Get(), dev_b.Get(), dev_c.Get());
   tryKernelLaunch();
   tryKernelSynchronize();
   auto stop = std::chrono::steady_clock::now();
   auto time = chrono::duration_cast<chrono::microseconds>(stop - start).count() * 1e-6;

   ////real time = timer.Ellapsed();
   //start = std::chrono::steady_clock::now();
   //AddMatricesShared<<<gridDim, blockDim>>>(dev_a.Get(), dev_b.Get(), dev_c.Get());
   //tryKernelLaunch();
   //tryKernelSynchronize();
   //stop = std::chrono::steady_clock::now();
   //shared_time = chrono::duration_cast<chrono::microseconds>(stop - start).count() * 1e-6;

   //real shared_time = timer.Ellapsed() - time;

   //long long lx = pow(2,8);
   //long long ly = pow(2,8);
   //long long lz = pow(2,8);

   //dim3 blockDim(16, 16, 4);
  // dim3 gridDim(lx / blockDim.x, ly / blockDim.y, lz / blockDim.z);

   //auto start = std::chrono::steady_clock::now();
   //Dummy<<<gridDim, blockDim >>>();
   //tryKernelLaunch();
   //tryKernelSynchronize();
   //auto stop = std::chrono::steady_clock::now();
   //auto dummy_time = chrono::duration_cast<chrono::microseconds>(stop - start).count() * 1e-6;

   //real dummy_time = timer.Ellapsed() - shared_time;

   //cout << setw(30) << "adding " << time << " " << 1 / time << endl;
   //cout << setw(30) << "adding with shared memory " << shared_time << " " << 1 / shared_time << endl;
   double runs_per_second = 1.0 / time;
   double tera = pow(2, 30);
   const double calc_per_thread = 2 * 20000.0;
   double TFLOPS = runs_per_second * calc_per_thread / tera * MATRIX_LENGTH;
   cout << setw(30) << "dummy " << time << " " << TFLOPS << endl;

   /*ofstream fout("results/c.txt");

   fout << scientific;

   for(size_t i = 0; i < MATRIX_WIDTH; i++)
   {
      for(size_t j = 0; j < MATRIX_WIDTH; j++)
      {
         fout << setw(16) << c[i * MATRIX_WIDTH + j] << " ";
      }

      fout << endl;
   }

   fout.close();*/

   /*delete a;
   delete b;
   delete c;*/
}