#include "cuda_timer.h"

using namespace cuda_utilities;

CudaTimer::CudaTimer()
{
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
}

CudaTimer::~CudaTimer()
{
   cudaEventDestroy(start);
   cudaEventDestroy(stop);
}

void CudaTimer::Start()
{
   cudaEventRecord(start, 0);
}

float CudaTimer::Ellapsed()
{
   float ellapsed_time = 0;

   cudaEventRecord(stop, 0);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&ellapsed_time, start, stop);

   return ellapsed_time * 1e-3;
}