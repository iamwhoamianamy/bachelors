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

public:
   const int n = 10;
   const int harmonicLength = (n + 1) * (n + 1);
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
   void calcLocalMultipolesWithLayersOrMatrices(bool useGPU, bool useMatrices);

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

   void calcMultipolesAtLeaves(
      const std::vector<std::vector<OctreeNode*>>& layers);

   Matrix<OctreeNode*> separateByOrientation(
      const std::vector<OctreeNode*>& layer);

   std::vector<ComplexMatrix> calcRegularMatrices(
      const Matrix<OctreeNode*>& nodesByOrientation);

   ComplexMatrix matrixFromRegularHarmonic(
      const ComplexHarmonicSeries& regular);

   std::vector<ComplexMatrix> getComponentsOfExpansionsInOneOrientation(
      const std::vector<OctreeNode*>& nodesByOrientation);

   void accountContributions(
      const std::vector<OctreeNode*>& nodesByOrientation,
      const std::vector<ComplexMatrix>& contributions);

};