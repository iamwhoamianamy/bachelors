#include "mesh.h"
#include "quad_points.h"
#include "vector3.cuh"
#include "dev_ptr.cuh"
#include "laplace_solver_interface.h"

using namespace cuda_utilities;
using namespace triangle_quadratures;

class LaplaceSolverArrays : public LaplaceSolver
{
public:
   LaplaceSolverArrays();
   void PrepareData(vector<Vector3>& points, Mesh& mesh, QuadPoints& quadPoints);
   vector<double>& SolveCPU();
   void CopyToDevice();
   void SolveGPU();
   vector<double>& GetResultGPU();

private:
   int quadraturesCount = 0;
   int trianglesCount = 0;
   int pointsCount = 0;
   int quadPointsOrder = 0;

   vector<double> quadratures_X;
   vector<double> quadratures_Y;
   vector<double> quadratures_Z;

   vector<double> normals_X;
   vector<double> normals_Y;
   vector<double> normals_Z;

   vector<double> points_X;
   vector<double> points_Y;
   vector<double> points_Z;

   vector<double> weights;
   vector<double> areas;
   vector<double> result;

   DevPtr<double> dev_quadratures_X;
   DevPtr<double> dev_quadratures_Y;
   DevPtr<double> dev_quadratures_Z;

   DevPtr<double> dev_normals_X;
   DevPtr<double> dev_normals_Y;
   DevPtr<double> dev_normals_Z;

   DevPtr<double> dev_points_X;
   DevPtr<double> dev_points_Y;
   DevPtr<double> dev_points_Z;

   DevPtr<double> dev_weights;
   DevPtr<double> dev_areas;
   DevPtr<double> dev_result;
};