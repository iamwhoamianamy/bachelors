#pragma once
#include "mesh.h"
#include "basis_quadratures.h"
#include "vector3.cuh"
#include "dev_ptr.h"
#include "laplace_solver_abstract.h"

using namespace cuda_utilities;
using namespace triangle_quadratures;

struct QuadPoint
{
   Vector3 quad;
   Vector3 normal;
   float weight;
};

class LaplaceSolverStructs : public LaplaceSolver
{
public:
   LaplaceSolverStructs();
   void PrepareData(vector<Vector3>& points, Mesh& mesh, BasisQuadratures& basisQuads);
   vector<float>& SolveCPU();
   void CopyToDevice();
   void SolveGPU();
   vector<float>& GetResultGPU();

private:
   int quadraturesCount = 0;
   int trianglesCount = 0;
   int pointsCount = 0;
   int quadraturesOrder = 0;

   vector<QuadPoint> quadPoints;
   vector<Vector3> points;

   vector<float> results;

   DevPtr<QuadPoint> dev_quadPoints;
   DevPtr<Vector3> dev_points;

   DevPtr<float> dev_results;
};