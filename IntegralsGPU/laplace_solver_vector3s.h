#pragma once
#include "real.h"
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
   void PrepareData(const vector<Vector3>& points, const Mesh& mesh, const BasisQuadratures& basisQuads);
   vector<real>& SolveCPU();
   void CopyToDevice();
   void SolveGPU();
   vector<real>& GetResultGPU();
   ~LaplaceSolverVector3s();
   //AlgorythmGPU algorythmGPU;

private:
   int quadraturesCount = 0;
   int trianglesCount = 0;
   int pointsCount = 0;
   int quadraturesOrder = 0;

   vector<Vector3> quadPoints;
   vector<Vector3> normals;
   vector<Vector3> points;

   vector<real> weights;
   vector<real> areas;
   vector<real> results;

   DevPtr<Vector3> dev_quadPoints;
   DevPtr<Vector3> dev_normals;
   DevPtr<Vector3> dev_points;

   DevPtr<real> dev_weights;
   DevPtr<real> dev_areas;
   DevPtr<real> dev_results;
};