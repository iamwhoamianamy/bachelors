#pragma once
#include <iostream>
#include <fstream>
#include <iomanip>

#include "laplace_data.cuh"
#include "laplace_solver_kernels.cuh"
#include "cuda_helper.h"
#include "real.h"

using namespace laplace_data;
using namespace laplace_solver_kernels;
using namespace cuda_utilities;

void addMatricesTest();
