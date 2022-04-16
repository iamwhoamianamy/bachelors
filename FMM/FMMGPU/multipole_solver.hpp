#pragma once
#include <vector>
#include "vector3.cuh"
#include "tetrahedron.hpp"
#include "basis_quadratures.hpp"
#include "harmonics.hpp"
#include "quadrature.hpp"
#include "octree.hpp"
#include "typedefs.hpp"

enum class TranslationAlgorithm
{
   CPU = 0,
   GPUSimple,
   GPUBlockForHarmonic
};

class MultipoleSolver
{
private:
   std::vector<Quadrature>& _quadratures;
   OctreeNode* octreeRoot;
   bool _multipolesAreReady = false;
   bool _log = true;

public:
   const int harmonicOrder = 10;
   const int harmonicLength = (harmonicOrder + 1) * (harmonicOrder + 1);
   const real eps = 1e-6;
   const size_t octreeLeafCapacity;

   MultipoleSolver(std::vector<Quadrature>& quadratures,
                   size_t octreeLeafCapacity = 1000);

   void calcLocalMultipolesWithoutTranslation();
   void calcLocalMultipolesWithComplexTranslation();
   void calcLocalMultipolesWithRealTranslation();
   void calcLocalMultipolesWithLayers(bool useGPU);
   void calcLocalMultipolesWithMatrices(bool useGPU);
   Vector3 calcA(real current, const Vector3& point);
   Vector3 calcB(real current, const Vector3& point);
   ~MultipoleSolver();

private:
   void calcLocalMultipolesWithLayersOrMatrices(
      bool useGPU,
      bool useMatrices);

   void enumerateNodes(
      OctreeNode* node,
      std::vector<std::vector<OctreeNode*>>& layers, 
      size_t currentLayerId);

   void calcContributionsToHigherLevels(
      const std::vector<std::vector<OctreeNode*>>& layers,
      bool useGPU,
      bool useMatrices);

   void calcContributionsToHigherLevels(
      const std::vector<std::vector<OctreeNode*>>& layers,
      bool useGPU);

   std::vector<Vector3> calcContributionsToHigherLevel(
      const std::vector<OctreeNode*>& layer,
      bool useGPU);

   void calcContributionsToHigherLevelsWithMatrices(
      const std::vector<std::vector<OctreeNode*>>& layers,
      bool useGPU);

   void calcContributionsToHigherLevelsWithMatricesCPU(
      const std::vector<std::vector<OctreeNode*>>& layers);

   void calcContributionsToHigherLevelsWithMatricesGPU(
      const std::vector<std::vector<OctreeNode*>>& layers);

   void calcMultipolesAtLeaves(
      const std::vector<std::vector<OctreeNode*>>& layers);

   Matrix<OctreeNode*> separateNodesByOrientation(
      const std::vector<OctreeNode*>& layer);

   std::vector<ComplexMatrix> calcRegularMatricesForLayer(
      const Matrix<OctreeNode*>& nodesByOrientation);

   ComplexMatrix calcRegularMatricesForLayerAsVectors(
      const Matrix<OctreeNode*>& nodesByOrientation);

   ComplexMatrix formMatrixFromRegularHarmonics(
      const ComplexHarmonicSeries& regular);

   std::vector<Complex> formMatrixFromRegularHarmonicsAsVectors(
      const ComplexHarmonicSeries& regular);

   std::vector<ComplexMatrix> getExpansionsInOneOrientation(
      const std::vector<OctreeNode*>& nodesByOrientation);

   ComplexMatrix getExpansionsInOneOrientationAsVectors(
      const std::vector<OctreeNode*>& nodesByOrientation);

   void accountChildrenContributions(
      const std::vector<OctreeNode*>& nodesByOrientation,
      const std::vector<ComplexMatrix>& contributions);

   void accountChildrenContributions(
      const std::vector<OctreeNode*>& nodesByOrientation,
      const ComplexMatrix& contributions);

   void printMatrices(
      const std::vector<ComplexMatrix>& regularMatrices,
      const std::vector<ComplexMatrix>& expansionMatrices);
};