#pragma once
#include "real.h"
#include <chrono>

#include "flops_tests.cuh"
#include "cuda_timer.h"

void addMatricesTest()
{
   const int MATRIX_LENGTH = MATRIX_WIDTH * MATRIX_HEIGHT;

   real* a = new real[MATRIX_LENGTH];
   real* b = new real[MATRIX_LENGTH];
   real* c = new real[MATRIX_LENGTH];

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

   DevPtr<real> dev_a(a, MATRIX_LENGTH);
   DevPtr<real> dev_b(b, MATRIX_LENGTH);
   DevPtr<real> dev_c(c, MATRIX_LENGTH);

   cout << "Done copying to GPU!" << endl;

   dim3 gridDim(MATRIX_WIDTH / BLOCK_SIZE, MATRIX_HEIGHT / BLOCK_SIZE);
   dim3 blockDim(BLOCK_SIZE, BLOCK_SIZE);

   auto start = std::chrono::steady_clock::now();
   AddMatricesShared<<<gridDim, blockDim >>>(dev_a.Get(), dev_b.Get(), dev_c.Get());
   tryKernelLaunch();
   tryKernelSynchronize();
   auto stop = std::chrono::steady_clock::now();
   auto time = chrono::duration_cast<chrono::microseconds>(stop - start).count() * 1e-6;
   
   real runs_per_second = 1.0 / time;
   real tera = pow(2, 40);
   const real calc_per_thread = 2 * 200000.0;
   real TFLOPS = runs_per_second * calc_per_thread / tera * MATRIX_LENGTH;
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

   delete a;
   delete b;
   delete c;
}