#pragma once
#include"cuda_runtime.h"
#include <cstdlib>
#include <stdio.h>

cudaError_t HANDLE_ERROR(cudaError_t funct_return);

cudaError_t TryCudaMalloc(void** dev_ptr, const size_t size);

