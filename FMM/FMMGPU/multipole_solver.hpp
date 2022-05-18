#pragma once
#include <vector>
#include "vector3.cuh"
#include "basis_quadratures.hpp"
#include "harmonics.hpp"
#include "quadrature.hpp"
#include "quadrature_octree.hpp"
#include "calculation_point_octree.hpp"
#include "typedefs.hpp"

enum class M2MAlg
{
   NoTranslation = 0,
   ComplexTranslation,
   RealTranslation,
   Layers,
   Matrices,
};

enum class M2MDevice
{
   CPU = 0,
   GPU,
   Adaptive
};


enum class M2LAlg
{
   NoTranslation = 0,
   ComplexTranslation,
   Matrices,
};

class MultipoleSolver
{
private:
   std::vector<Quadrature>& _quadratures;
   std::vector<Vector3>& _points;
   QuadratureOctreeNode* _quadratureOctreeRoot;
   CalculationPointOctreeNode* _calculationPointOctreeRoot;
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
   const size_t octreeLeafCapacity;

   MultipoleSolver(
      std::vector<Quadrature>& quadratures,
      size_t octreeLeafCapacity = 1000);

   MultipoleSolver(
      std::vector<Vector3>& points,
      std::vector<Quadrature>& quadratures,
      size_t octreeLeafCapacity = 1000);

   void calcMultipoleExpansionsAtLeaves();

   size_t getOctreeNodeCount() const;

   void calclMultipoleExpansions(
      M2MAlg algorithm,
      M2MDevice device = M2MDevice::CPU);

   void calclLocalMultipoleExpansions(
      M2LAlg algorithm,
      M2MDevice device = M2MDevice::CPU);

   Vector3 calcA(real current, const Vector3& point);
   Vector3 calcB(real current, const Vector3& point);
   ~MultipoleSolver();

private:
   void calcMultipoleExpansionsWithoutTranslation();
   void calcMultipoleExpansionsWithComplexTranslation();
   void calcMultipoleExpansionsWithRealTranslation();

   void calcMultipoleExpanstionsWithLayersOrMatrices(
      M2MDevice device,
      bool useMatrices);

   void enumerateNodes(
      QuadratureOctreeNode* node,
      std::vector<std::vector<QuadratureOctreeNode*>>& layers, 
      size_t currentLayerId);

   void calcContributionsToHigherLayers(
      const std::vector<std::vector<QuadratureOctreeNode*>>& layers,
      M2MDevice device,
      bool useMatrices);

   void calcContributionsToHigherLayers(
      const std::vector<std::vector<QuadratureOctreeNode*>>& layers,
      M2MDevice device);

   std::vector<Vector3> calcContributionsToHigherLayer(
      const std::vector<QuadratureOctreeNode*>& layer,
      M2MDevice device);

   void calcContributionsToHigherLevelsWithMatrices(
      const std::vector<std::vector<QuadratureOctreeNode*>>& layers,
      M2MDevice device);

   void calcMultipoleExpansionsAtLeaves(
      const std::vector<std::vector<QuadratureOctreeNode*>>& layers);

   Matrix<QuadratureOctreeNode*> separateNodesByOrientation(
      const std::vector<QuadratureOctreeNode*>& layer);

   std::vector<ComplexMatrix> calcRegularMatricesForLayer(
      const Matrix<QuadratureOctreeNode*>& nodesByOrientation);

   RealMatrix calcRegularMatricesForLayerAsVectors(
      const Matrix<QuadratureOctreeNode*>& nodesByOrientation);

   ComplexMatrix formMatrixFromRegularHarmonics(
      const ComplexHarmonicSeries& regular);

   std::vector<Complex> formMatrixFromRegularHarmonicsAsVectors(
      const ComplexHarmonicSeries& regular);

   RealMatrix getExpansionsInOneOrientationAsVectors(
      const std::vector<QuadratureOctreeNode*>& nodesByOrientation);

   void accountChildrenContributions(
      const std::vector<QuadratureOctreeNode*>& nodesByOrientation,
      const RealMatrix& contributions);

   void printMatrices(
      const std::vector<ComplexMatrix>& regularMatrices,
      const std::vector<ComplexMatrix>& expansionMatrices);

   void initTrees(
      std::vector<Vector3>& points,
      std::vector<Quadrature>& quadratures,
      size_t octreeLeafCapacity);

   void initTransitionMatrices();

   void calcLocalMultipoleExpansionsWithoutTranslation();
   void calcLocalMultipoleExpansionsWithComplexTranslation();
};