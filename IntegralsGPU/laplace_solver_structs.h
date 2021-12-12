#pragma once
#include "real.h"
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
   real weight;
};

class LaplaceSolverStructs : public LaplaceSolver
{
public:
   LaplaceSolverStructs();
   void PrepareData(const vector<Vector3>& points,
                    const Mesh& mesh,
                    const BasisQuadratures& basisQuads);
   vector<real>& SolveCPU();
   void CopyToDevice();
   void SolveGPU();
   vector<real>& GetResultGPU();
   size_t PointsCountPadded() const;
   ~LaplaceSolverStructs();

private:
   size_t quadraturesCount = 0;
   size_t pointsCount = 0;
   size_t matrixWidth = 0;

   vector<QuadPoint> quadPoints;
   vector<Vector3> points;
   vector<real> results;

   DevPtr<QuadPoint> dev_quadPoints;
   DevPtr<Vector3> dev_points;
   DevPtr<real> dev_results;
};