#pragma once
#include "mesh.h"
#include "basis_quadratures.h"
#include "vector3.cuh"

using namespace triangle_quadratures;

enum class AlgorythmGPU
{
   Reduction,
   Blocks
};

class LaplaceSolver
{
public:
   virtual void PrepareData(vector<Vector3>& points, Mesh& mesh, BasisQuadratures& basisQuads) = 0;
   virtual vector<float>& SolveCPU() = 0;
   virtual void CopyToDevice() = 0;
   virtual void SolveGPU() = 0;
   virtual vector<float>& GetResultGPU() = 0;

   AlgorythmGPU algorythmGPU;
};