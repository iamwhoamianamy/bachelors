#pragma once
#include "mesh.h"
#include "quad_points.h"
#include "vector3.cuh"

using namespace triangle_quadratures;

class LaplaceSolver
{
public:
   virtual void PrepareData(vector<Vector3>& points, Mesh& mesh, QuadPoints& quadPoints) = 0;
   virtual vector<double>& SolveCPU() = 0;
   virtual void CopyToDevice() = 0;
   virtual void SolveGPU() = 0;
   virtual vector<double>& GetResultGPU() = 0;
};