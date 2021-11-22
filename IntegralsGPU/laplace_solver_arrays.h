#pragma once
#include "real.h"
#include "mesh.h"
#include "basis_quadratures.h"
#include "vector3.cuh"
#include "dev_ptr.h"
#include "laplace_solver_abstract.h"

using namespace cuda_utilities;
using namespace triangle_quadratures;

class LaplaceSolverArrays : public LaplaceSolver
{
public:
   LaplaceSolverArrays();
   void PrepareData(const vector<Vector3>& points, const Mesh& mesh, const BasisQuadratures& basisQuads);
   vector<real>& SolveCPU();
   void CopyToDevice();
   void SolveGPU();
   vector<real>& GetResultGPU();
   ~LaplaceSolverArrays();

private:
   int quadraturesCount = 0;
   int trianglesCount = 0;
   int pointsCount = 0;
   int quadraturesOrder = 0;

   vector<real> quadratures_X;
   vector<real> quadratures_Y;
   vector<real> quadratures_Z;

   vector<real> normals_X;
   vector<real> normals_Y;
   vector<real> normals_Z;

   vector<real> points_X;
   vector<real> points_Y;
   vector<real> points_Z;

   vector<real> weights;
   vector<real> areas;
   vector<real> results;

   DevPtr<real> dev_quadratures_X;
   DevPtr<real> dev_quadratures_Y;
   DevPtr<real> dev_quadratures_Z;

   DevPtr<real> dev_normals_X;
   DevPtr<real> dev_normals_Y;
   DevPtr<real> dev_normals_Z;

   DevPtr<real> dev_points_X;
   DevPtr<real> dev_points_Y;
   DevPtr<real> dev_points_Z;

   DevPtr<real> dev_weights;
   DevPtr<real> dev_areas;
   DevPtr<real> dev_results;
};