#pragma once
#include <vector>
#include "vector3.cuh"
#include "tetrahedron.hpp"
#include "basis_quadratures.hpp"
#include "harmonics.hpp"
#include "quadrature.hpp"
#include "octree.hpp"
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

class MultipoleSolver
{
private:
   std::vector<Quadrature>& _quadratures;
   OctreeNode* octreeRoot;
   bool _multipolesAreReady = false;

public:
   bool log = true;
   float adaptiveBorder = 16e3;
   const int harmonicOrder = 10;
   const int harmonicLength = (harmonicOrder + 1) * (harmonicOrder + 1);
   const real eps = 1e-6;
   const size_t octreeLeafCapacity;

   MultipoleSolver(std::vector<Quadrature>& quadratures,
                   size_t octreeLeafCapacity = 1000);

   void calcLocalMultipoles(
      M2MAlg algorithm,
      M2MDevice device = M2MDevice::CPU);

   Vector3 calcA(real current, const Vector3& point);
   Vector3 calcB(real current, const Vector3& point);
   ~MultipoleSolver();

private:
   void calcLocalMultipolesWithoutTranslation();
   void calcLocalMultipolesWithComplexTranslation();
   void calcLocalMultipolesWithRealTranslation();

   void calcLocalMultipolesWithLayersOrMatrices(
      M2MDevice device,
      bool useMatrices);

   void enumerateNodes(
      OctreeNode* node,
      std::vector<std::vector<OctreeNode*>>& layers, 
      size_t currentLayerId);

   void calcContributionsToHigherLayers(
      const std::vector<std::vector<OctreeNode*>>& layers,
      M2MDevice device,
      bool useMatrices);

   void calcContributionsToHigherLayers(
      const std::vector<std::vector<OctreeNode*>>& layers,
      M2MDevice device);

   std::vector<Vector3> calcContributionsToHigherLayer(
      const std::vector<OctreeNode*>& layer,
      M2MDevice device);

   void calcContributionsToHigherLevelsWithMatrices(
      const std::vector<std::vector<OctreeNode*>>& layers,
      M2MDevice device);

   void calcMultipolesAtLeaves(
      const std::vector<std::vector<OctreeNode*>>& layers);

   Matrix<OctreeNode*> separateNodesByOrientation(
      const std::vector<OctreeNode*>& layer);

   std::vector<ComplexMatrix> calcRegularMatricesForLayer(
      const Matrix<OctreeNode*>& nodesByOrientation);

   ComplexMatrix calcRegularMatricesForLayerAsVectors(
      const Matrix<OctreeNode*>& nodesByOrientation,
      size_t padding = 0);

   ComplexMatrix formMatrixFromRegularHarmonics(
      const ComplexHarmonicSeries& regular);

   std::vector<Complex> formMatrixFromRegularHarmonicsAsVectors(
      const ComplexHarmonicSeries& regular,
      size_t padding = 0);

   std::vector<ComplexMatrix> getExpansionsInOneOrientation(
      const std::vector<OctreeNode*>& nodesByOrientation);

   ComplexMatrix getExpansionsInOneOrientationAsVectors(
      const std::vector<OctreeNode*>& nodesByOrientation,
      size_t padding = 0);

   void accountChildrenContributions(
      const std::vector<OctreeNode*>& nodesByOrientation,
      const std::vector<ComplexMatrix>& contributions);

   void accountChildrenContributions(
      const std::vector<OctreeNode*>& nodesByOrientation,
      const ComplexMatrix& contributions,
      size_t padding = 0);

   void printMatrices(
      const std::vector<ComplexMatrix>& regularMatrices,
      const std::vector<ComplexMatrix>& expansionMatrices);
};