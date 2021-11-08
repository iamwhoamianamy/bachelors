#include "mesh.h"
#include "quad_points.h"
#include "vector3.cuh"
#include "dev_ptr.cuh"
#include "laplace_solver_interface.h"

using namespace cuda_utilities;
using namespace triangle_quadratures;

class LaplaceSolverVector3s : public LaplaceSolver
{
public:
   LaplaceSolverVector3s();
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

   vector<Vector3> quadratures;
   vector<Vector3> normals;
   vector<Vector3> points;

   vector<double> weights;
   vector<double> areas;
   vector<double> result;

   DevPtr<Vector3> dev_quadratures;
   DevPtr<Vector3> dev_normals;
   DevPtr<Vector3> dev_points;

   DevPtr<double> dev_weights;
   DevPtr<double> dev_areas;
   DevPtr<double> dev_result;
};