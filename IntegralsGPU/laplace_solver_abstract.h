#pragma once
#include "real.h"
#include "mesh.h"
#include "basis_quadratures.h"
#include "vector3.cuh"

using namespace triangle_quadratures;

enum class AlgorythmGPU
{
   Reduction,
   Blocks,
   Grid
};

class LaplaceSolver
{
public:
   virtual void PrepareData(const vector<Vector3>& points,
                            const Mesh& mesh,
                            const BasisQuadratures& basisQuads) = 0;
   virtual vector<real>& SolveCPU() = 0;
   virtual void CopyToDevice() = 0;
   virtual void SolveGPU() = 0;
   virtual vector<real>& GetResultGPU() = 0;
   virtual ~LaplaceSolver() = 0;
};
