#pragma once
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

   vector<float> quadratures_X;
   vector<float> quadratures_Y;
   vector<float> quadratures_Z;

   vector<float> normals_X;
   vector<float> normals_Y;
   vector<float> normals_Z;

   vector<float> points_X;
   vector<float> points_Y;
   vector<float> points_Z;

   vector<float> weights;
   vector<float> areas;
   vector<float> results;

   DevPtr<float> dev_quadratures_X;
   DevPtr<float> dev_quadratures_Y;
   DevPtr<float> dev_quadratures_Z;

   DevPtr<float> dev_normals_X;
   DevPtr<float> dev_normals_Y;
   DevPtr<float> dev_normals_Z;

   DevPtr<float> dev_points_X;
   DevPtr<float> dev_points_Y;
   DevPtr<float> dev_points_Z;

   DevPtr<float> dev_weights;
   DevPtr<float> dev_areas;
   DevPtr<float> dev_results;
};