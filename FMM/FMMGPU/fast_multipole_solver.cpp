#include "fast_multipole_solver.hpp"

FastMultipoleSolver::FastMultipoleSolver(
   std::vector<Quadrature>& quadratures,
   std::vector<Vector3>& points,
   size_t quadratureOctreeLeafCapacity,
   size_t calculationPointOctreeLeafCapacity) :
   MultipoleSolver(quadratures, quadratureOctreeLeafCapacity),
   _points(points),
   calculationPointOctreeLeafCapacity(calculationPointOctreeLeafCapacity)
{
   initTrees();
}

FastMultipoleSolver::~FastMultipoleSolver()
{
   delete _calculationPointOctreeRoot;
}

void FastMultipoleSolver::initTrees()
{
   MultipoleSolver::initTrees();

   _calculationPointOctreeRoot = new CalculationPointOctreeNode(
      //Box(Vector3(0, 0, 0), Vector3(3, 3, 3)), calculationPointOctreeLeafCapacity);
      Box(Vector3(10, 5, 8), Vector3(1, 1, 1)), calculationPointOctreeLeafCapacity);
   _calculationPointOctreeRoot->insert(_points);
}

void FastMultipoleSolver::calcLocalMultipoleExpansionsWithComplexTranslation() const
{
   auto nodesToVisit = _calculationPointOctreeRoot->getAllNodesAsSet();

   _quadratureOctreeRoot->translateMultipoleExpansionsToLocal(
      _calculationPointOctreeRoot,
      nodesToVisit);

   _calculationPointOctreeRoot->propagateLocalExpansions();
}
void FastMultipoleSolver::calclLocalMultipoleExpansions(
   M2LAlg algorithm, 
   M2MDevice device)
{
   if(_multipolesAreReady)
   {
      if(!_localMultipolesAreInitialized)
      {
         _calculationPointOctreeRoot->initAllLocalExpansions(harmonicOrder);
         _localMultipolesAreInitialized = true;
      }

      switch(algorithm)
      {
         case M2LAlg::ComplexTranslation:
         {
            calcLocalMultipoleExpansionsWithComplexTranslation();
            break;
         }
         default:
         {
            throw std::exception("Not implemented!");
         }
      }
   }
   else
   {
      throw std::exception("Multipoles are not ready!");
   }
}

Vector3 FastMultipoleSolver::calcA(real current, const Vector3& point)
{
   return MultipoleSolver::calcA(current, point);
}

Vector3 FastMultipoleSolver::calcB(real current, const Vector3& point)
{
   return MultipoleSolver::calcB(current, point);
}

std::vector<std::pair<Vector3, Vector3>> FastMultipoleSolver::calcA(real current)
{
   auto res = _calculationPointOctreeRoot->calcA(_points.size());

   for (auto &re : res)
   {
      re.second = re.second / (4.0 * math::PI) * current;
   }

   return res;
}

std::vector<Vector3> FastMultipoleSolver::calcB(real current)
{
   std::vector<Vector3> res(_points.size());



   return res;
}
