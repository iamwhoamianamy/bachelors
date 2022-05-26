#pragma once
#include <vector>
#include "vector3.cuh"
#include "basis_quadratures.hpp"
#include "harmonics.hpp"
#include "quadrature.hpp"
#include "quadrature_octree.hpp"
#include "typedefs.hpp"
#include "multipole_solver_enums.hpp"

class MultipoleSolver
{
protected:
   std::vector<Quadrature>& _quadratures;
   QuadratureOctreeNode* _quadratureOctreeRoot = nullptr;
   bool _multipolesAreReady = false;
   bool _multipolesAtLeavesAreReady = false;
   std::vector<Complex> _realToComplexMatrix;
   std::vector<Complex> _complexToRealMatrix;

public:
   bool log = true;
   float adaptiveBorder = 16e3;
   const int harmonicOrder = 10;
   const int harmonicLength = (harmonicOrder + 1) * (harmonicOrder + 1);
   const real eps = 1e-6;
   const size_t quadratureOctreeLeafCapacity;

   MultipoleSolver(
      std::vector<Quadrature>& quadratures,
      size_t quadratureOctreeLeafCapacity = 1000);

   virtual void calcMultipoleExpansionsAtLeaves();

   virtual size_t getOctreeNodeCount() const;

   virtual void calclMultipoleExpansions(
      M2MAlg algorithm,
      Device device = Device::CPU);
   
   virtual Vector3 calcA(real current, const Vector3& point);
   virtual Vector3 calcB(real current, const Vector3& point);
   virtual ~MultipoleSolver();

protected:
   virtual void calcMultipoleExpansionsWithoutTranslation();
   virtual void calcMultipoleExpansionsWithComplexTranslation();
   virtual void calcMultipoleExpansionsWithRealTranslation();

   virtual void calcMultipoleExpanstionsWithLayersOrMatrices(
      Device device,
      bool useMatrices);

   virtual void enumerateNodes(
      QuadratureOctreeNode* node,
      std::vector<std::vector<QuadratureOctreeNode*>>& layers, 
      size_t currentLayerId);

   virtual void calcContributionsToHigherLayers(
      const std::vector<std::vector<QuadratureOctreeNode*>>& layers,
      Device device,
      bool useMatrices);

   virtual void calcContributionsToHigherLayers(
      const std::vector<std::vector<QuadratureOctreeNode*>>& layers,
      Device device);

   virtual std::vector<Vector3> calcContributionsToHigherLayer(
      const std::vector<QuadratureOctreeNode*>& layer,
      Device device);

   virtual void calcContributionsToHigherLevelsWithMatrices(
      const std::vector<std::vector<QuadratureOctreeNode*>>& layers,
      Device device);

   virtual void calcMultipoleExpansionsAtLeaves(
      const std::vector<std::vector<QuadratureOctreeNode*>>& layers);

   virtual Matrix<QuadratureOctreeNode*> separateNodesByOrientation(
      const std::vector<QuadratureOctreeNode*>& layer);

   virtual std::vector<ComplexMatrix> calcRegularMatricesForLayer(
      const Matrix<QuadratureOctreeNode*>& nodesByOrientation);

   virtual RealMatrix calcRegularMatricesForLayerAsVectors(
      const Matrix<QuadratureOctreeNode*>& nodesByOrientation);

   virtual ComplexMatrix formMatrixFromRegularHarmonics(
      const ComplexHarmonicSeries& regular);

   virtual std::vector<Complex> formMatrixFromRegularHarmonicsAsVectors(
      const ComplexHarmonicSeries& regular);

   virtual RealMatrix getExpansionsInOneOrientationAsVectors(
      const std::vector<QuadratureOctreeNode*>& nodesByOrientation);

   virtual void accountChildrenContributions(
      const std::vector<QuadratureOctreeNode*>& nodesByOrientation,
      const RealMatrix& contributions);

   virtual void printMatrices(
      const std::vector<ComplexMatrix>& regularMatrices,
      const std::vector<ComplexMatrix>& expansionMatrices);

   virtual void initTrees();
   virtual void initTransitionMatrices();
};