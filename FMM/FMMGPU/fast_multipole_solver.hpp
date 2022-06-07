#pragma once
#include <map>
#include <vector>
#include "multipole_solver.hpp"
#include "vector3.cuh"
#include "basis_quadratures.hpp"
#include "quadrature.hpp"
#include "calculation_point_octree.hpp"
#include "multipole_solver_enums.hpp"

typedef std::map<CalculationPointOctreeNode*, std::set<QuadratureOctreeNode*>> InteractionMap;

class FastMultipoleSolver : public MultipoleSolver
{
private:
   std::vector<Vector3>& _points;
   CalculationPointOctreeNode* _calculationPointOctreeRoot = nullptr;
   bool _localMultipolesAreInitialized = false;
   InteractionMap _farInteractionMap;
   InteractionMap _closeInteractionMap;

public:
   const size_t calculationPointOctreeLeafCapacity;

   FastMultipoleSolver(
      std::vector<Quadrature>& quadratures,
      std::vector<Vector3>& points,
      size_t quadratureOctreeLeafCapacity = 1000,
      size_t calculationPointOctreeLeafCapacity = 100);
   
   void calcLocalMultipoleExpansions(
      M2LAlg algorithm,
      Device device = Device::CPU);
   
   virtual Vector3 calcA(real current, const Vector3& point) override;
   virtual Vector3 calcB(real current, const Vector3& point) override;
   virtual std::vector<std::pair<Vector3, Vector3>> calcA(real current);
   virtual std::vector<std::pair<Vector3, Vector3>> calcB(real current);

   virtual ~FastMultipoleSolver() override;

private:
   void initTrees() override;
   void propagateLocalExpansionsWithComplexTranslation();
   void propagateLocalExpansionsWithMatrices(Device device);

   void formInteractionMaps(
      QuadratureOctreeNode* quadratureTreeNode,
      CalculationPointOctreeNode* calcPointTreeNode);

   bool checkIfFarEnough(
      const QuadratureOctreeNode* quadratureTreeNode,
      const CalculationPointOctreeNode* calcPointTreeNode) const;

   void accountFarInteractions();

   void calcContributionsToLowerLayers(
      const std::vector<std::vector<CalculationPointOctreeNode*>>& layers,
      Device device);

   void accountContributionsToChildren(
      const std::vector<CalculationPointOctreeNode*>& nodesByOrientation,
      const RealMatrix& contributions);

   RealMatrix getExpansionsInOneOrientationAsVectors(
      const std::vector<CalculationPointOctreeNode*>& nodesByOrientation) const;

   virtual void enumerateNodes(
      QuadratureOctreeNode* node,
      std::vector<std::vector<QuadratureOctreeNode*>>& layers,
      size_t currentLayerId) override;

   virtual void enumerateNodes(
      CalculationPointOctreeNode* node,
      std::vector<std::vector<CalculationPointOctreeNode*>>& layers,
      size_t currentLayerId);

   virtual Matrix<QuadratureOctreeNode*> separateNodesByOrientation(
      const std::vector<QuadratureOctreeNode*>& layer) override;

   Matrix<CalculationPointOctreeNode*> separateNodesByOrientation(
      const std::vector<CalculationPointOctreeNode*>& layer);

   RealMatrix calcRegularMatricesForL2LAsVectors(
      const Matrix<CalculationPointOctreeNode*>& nodesByOrientation);

   std::vector<Complex> formMatrixFromRegularHarmonicsForL2LAsVectors(
      const ComplexHarmonicSeries& regular) const;
};