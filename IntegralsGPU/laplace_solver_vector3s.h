#pragma once
#include "mesh.h"
#include "basis_quadratures.h"
#include "vector3.cuh"
#include "dev_ptr.h"
#include "laplace_solver_abstract.h"

using namespace cuda_utilities;
using namespace triangle_quadratures;

class LaplaceSolverVector3s : public LaplaceSolver
{
public:
   LaplaceSolverVector3s();
   void PrepareData(vector<Vector3>& points, Mesh& mesh, BasisQuadratures& basisQuads);
   vector<float>& SolveCPU();
   void CopyToDevice();
   void SolveGPU();
   vector<float>& GetResultGPU();
   //AlgorythmGPU algorythmGPU;

private:
   int quadraturesCount = 0;
   int trianglesCount = 0;
   int pointsCount = 0;
   int quadraturesOrder = 0;

   vector<Vector3> quadPoints;
   vector<Vector3> normals;
   vector<Vector3> points;

   vector<float> weights;
   vector<float> areas;
   vector<float> results;

   DevPtr<Vector3> dev_quadPoints;
   DevPtr<Vector3> dev_normals;
   DevPtr<Vector3> dev_points;

   DevPtr<float> dev_weights;
   DevPtr<float> dev_areas;
   DevPtr<float> dev_results;
};