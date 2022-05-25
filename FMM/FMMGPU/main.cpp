#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>

#include "cblas.h"
#include "vector3.cuh"
#include "math.hpp"
#include "integration.hpp"
#include "harmonics.hpp"
#include "multipole_solver.hpp"
#include "fast_multipole_solver.hpp"
#include "translation_algorithms.hpp"
#include "testing_helpers.hpp"
#include "multipole_translator.hpp"

using namespace math;
using namespace test;
constexpr real current = 5;

std::vector<std::pair<Vector3, Vector3>> readTelmaResults(const std::string& filename)
{
   std::vector<std::pair<Vector3, Vector3>> res;
   std::ifstream fin(filename);
   size_t pointsCount;
   std::string _;
   fin >> _ >> _ >> _;
   fin >> pointsCount;
   res.resize(pointsCount);

   for(size_t i = 0; i < pointsCount; i++)
   {
      real px, py, pz, hx, hy, hz;
      
      fin >> _ >> px >> py >> pz >> hx >> hy >> hz;

      res[i] = std::pair<Vector3, Vector3>(Vector3(px, py, pz), Vector3(hx, hy, hz));
   }

   return res;
}

void comparisonToTelmaIntegrals()
{
   Torus torus = createTorus();
   BasisQuadratures bq = readBasisQuadratures();
   auto telmaResults = readTelmaResults("results/telma_results.txt");
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

   MultipoleSolver multipoleSolver(quadratures);
   multipoleSolver.calcMultipoleExpansionsAtLeaves();
   multipoleSolver.log = false;

   multipoleSolver.calclMultipoleExpansions(M2MAlg::NoTranslation);
   multipoleSolver.calclMultipoleExpansions(M2MAlg::ComplexTranslation);
   multipoleSolver.calclMultipoleExpansions(M2MAlg::RealTranslation);
   multipoleSolver.calclMultipoleExpansions(M2MAlg::Layers, M2MDevice::CPU);
   multipoleSolver.calclMultipoleExpansions(M2MAlg::Layers, M2MDevice::GPU);
   multipoleSolver.calclMultipoleExpansions(M2MAlg::Matrices, M2MDevice::CPU);
   multipoleSolver.calclMultipoleExpansions(M2MAlg::Matrices, M2MDevice::GPU);

   real sumErrorIntegration = 0;
   real sumErrorNoTranslation = 0;
   real sumErrorComplexTranslation = 0;
   real sumErrorRealTranslation = 0;
   real sumErrorLayersTranslationCPU = 0;
   real sumErrorLayersTranslationGPU = 0;
   real sumErrorMatricesTranslationCPU = 0;
   real sumErrorMatricesTranslationGPU = 0;

   size_t n = telmaResults.size();

   for(size_t i = 0; i < n; i++)
   {
      auto point = telmaResults[i].first;
      auto telmaB = telmaResults[i].second * math::MU0;

      Vector3 integration = math::calcBioSavartLaplace(current, point, quadratures);
      Vector3 noTranslation = multipoleSolver.calcB(current, point);
      Vector3 complexTranslation = multipoleSolver.calcB(current, point);
      Vector3 realTranslation = multipoleSolver.calcB(current, point);
      Vector3 layersTranslationCPU = multipoleSolver.calcB(current, point);
      Vector3 layersTranslationGPU = multipoleSolver.calcB(current, point);
      Vector3 matricesTranslationCPU = multipoleSolver.calcB(current, point);
      Vector3 matricesTranslationGPU = multipoleSolver.calcB(current, point);

      real errorIntegration = 100 * (integration - telmaB).length() / telmaB.length();
      real errorNoTranslation = 100 * (noTranslation - telmaB).length() / telmaB.length();
      real errorComplexTranslation = 100 * (complexTranslation - telmaB).length() / telmaB.length();
      real errorRealTranslation = 100 * (realTranslation - telmaB).length() / telmaB.length();
      real errorLayersTranslationCPU = 100 * (layersTranslationCPU - telmaB).length() / telmaB.length();
      real errorLayersTranslationGPU = 100 * (layersTranslationGPU - telmaB).length() / telmaB.length();
      real errorMatricesTranslationCPU = 100 * (matricesTranslationCPU - telmaB).length() / telmaB.length();
      real errorMatricesTranslationGPU = 100 * (matricesTranslationGPU - telmaB).length() / telmaB.length();

      std::cout << std::fixed << std::setw(3) << i << " ";
      point.printWithWidth(std::cout, 6);
      std::cout << std::setw(16) << errorIntegration;
      std::cout << std::setw(16) << errorNoTranslation;
      std::cout << std::setw(16) << errorComplexTranslation;
      std::cout << std::setw(16) << errorRealTranslation;
      std::cout << std::setw(16) << errorLayersTranslationCPU;
      std::cout << std::setw(16) << errorLayersTranslationGPU;
      std::cout << std::setw(16) << errorMatricesTranslationCPU;
      std::cout << std::setw(16) << errorMatricesTranslationGPU;
      std::cout << std::endl;

      sumErrorNoTranslation += errorNoTranslation;
      sumErrorComplexTranslation += errorComplexTranslation;
      sumErrorRealTranslation += errorRealTranslation;
      sumErrorLayersTranslationCPU += errorLayersTranslationCPU;
      sumErrorLayersTranslationGPU += errorLayersTranslationGPU;
      sumErrorMatricesTranslationCPU += errorMatricesTranslationCPU;
      sumErrorMatricesTranslationGPU += errorMatricesTranslationGPU;
   }

   std::cout << sumErrorNoTranslation / n << std::endl;
   std::cout << sumErrorComplexTranslation / n << std::endl;
   std::cout << sumErrorRealTranslation / n << std::endl;
   std::cout << sumErrorLayersTranslationCPU / n << std::endl;
   std::cout << sumErrorLayersTranslationGPU / n << std::endl;
   std::cout << sumErrorMatricesTranslationCPU / n << std::endl;
   std::cout << sumErrorMatricesTranslationGPU / n << std::endl;
}

void comparisonBetweenMethodsOnPrecision()
{
   Torus torus = createTorus();
   BasisQuadratures bq = readBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);
   MultipoleSolver multipoleSolver(quadratures, 10);
   multipoleSolver.log = false;
   multipoleSolver.calcMultipoleExpansionsAtLeaves();

   Vector3 point(10, 5, 8);

   Vector3 byIntegration = math::calcBioSavartLaplace(current, point, quadratures);
   
   multipoleSolver.calclMultipoleExpansions(M2MAlg::NoTranslation);
   Vector3 byMultipolesWithoutTranslation = multipoleSolver.calcB(current, point);

   multipoleSolver.calclMultipoleExpansions(M2MAlg::ComplexTranslation);
   Vector3 byMultipolesWithComplexTranslation = multipoleSolver.calcB(current, point);

   multipoleSolver.calclMultipoleExpansions(M2MAlg::RealTranslation);
   Vector3 byMultipolesWithRealTranslation = multipoleSolver.calcB(current, point);

   multipoleSolver.calclMultipoleExpansions(M2MAlg::Layers);
   Vector3 byMultipolesWithLayersCPU = multipoleSolver.calcB(current, point);

   multipoleSolver.calclMultipoleExpansions(M2MAlg::Layers, M2MDevice::GPU);
   Vector3 byMultipolesWithLayersGPU = multipoleSolver.calcB(current, point);

   multipoleSolver.calclMultipoleExpansions(M2MAlg::Matrices, M2MDevice::CPU);
   Vector3 byMultipolesWithMatricesCPU = multipoleSolver.calcB(current, point);

   multipoleSolver.calclMultipoleExpansions(M2MAlg::Matrices, M2MDevice::GPU);
   Vector3 byMultipolesWithMatricesGPU = multipoleSolver.calcB(current, point);

   std::cout << std::setw(20) << "point ";
   point.printWithWidth(std::cout, 6);
   std::cout << std::scientific << std::endl;
   std::cout << std::setw(40) << "integration " << byIntegration << std::endl;
   std::cout << std::setw(40) << "multipoles w/t translation " << byMultipolesWithoutTranslation << std::endl;
   std::cout << std::setw(40) << "multipoles with c translation " << byMultipolesWithComplexTranslation << std::endl;
   std::cout << std::setw(40) << "multipoles with r translation " << byMultipolesWithRealTranslation << std::endl;
   std::cout << std::setw(40) << "multipoles with layers CPU" << byMultipolesWithLayersCPU << std::endl;
   std::cout << std::setw(40) << "multipoles with layers GPU" << byMultipolesWithLayersGPU << std::endl;
   std::cout << std::setw(40) << "multipoles with matrices CPU" << byMultipolesWithMatricesCPU << std::endl;
   std::cout << std::setw(40) << "multipoles with matrices GPU" << byMultipolesWithMatricesGPU << std::endl;
}

void octreeFormingTime()
{
   constexpr double torusRadius = 2;
   constexpr double torusSectionWidth = 0.2;

   size_t w = 15;

   for(size_t i = 1; i < 16; i++)
   {
      Torus torus(torusRadius, torusSectionWidth, 10 * i, 16, 16);
      auto bq = readBasisQuadratures();
      auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

      auto start = std::chrono::steady_clock::now();
      MultipoleSolver multipoleSolver(quadratures, 1000);
      auto stop = std::chrono::steady_clock::now();
      double timeForOctreeForming = getTime(start, stop);

      std::cout << quadratures.size() << " ";
      std::cout << timeForOctreeForming << std::endl;
   }
}

void calculationTimeForMultipolesInLeaves()
{
   constexpr double torusRadius = 2;
   constexpr double torusSectionWidth = 0.2;

   size_t w = 15;

   for(size_t i = 1; i < 15; i++)
   {
      Torus torus(torusRadius, torusSectionWidth, 10 * i, 16, 16);
      auto bq = readBasisQuadratures();
      auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

      MultipoleSolver multipoleSolver(quadratures, 16);

      auto start = std::chrono::steady_clock::now();
      multipoleSolver.calcMultipoleExpansionsAtLeaves();
      auto stop = std::chrono::steady_clock::now();
      double timeForMultipolesInLeaves = getTime(start, stop);

      std::cout << multipoleSolver.getOctreeNodeCount() << " ";
      std::cout << timeForMultipolesInLeaves << std::endl;
   }
}

void calculationTimeForLocalMultipolesByLeafCapacity()
{
   auto torus = createTorus();
   auto bq = readBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

   size_t w = 15;

   std::cout << std::setw(w) << "leaf capacity";
   std::cout << std::setw(w) << "w/t";
   std::cout << std::setw(w) << "Complex";
   std::cout << std::setw(w) << "real";
   std::cout << std::setw(w) << "layersCPU";
   std::cout << std::setw(w) << "layersGPU";
   std::cout << std::setw(w) << "matricesCPU";
   std::cout << std::setw(w) << "matricesGPU";
   //std::cout << std::setw(w) << "matricesAda";
   std::cout << std::endl;

   std::cout << std::fixed;

   for(size_t i = 3; i < 15; i++)
   {
      int octreeLeafCapacity = pow(2, i);
      MultipoleSolver multipoleSolver(quadratures, octreeLeafCapacity);
      multipoleSolver.calcMultipoleExpansionsAtLeaves();
      multipoleSolver.log = false;

      auto start = std::chrono::steady_clock::now();
      multipoleSolver.calclMultipoleExpansions(M2MAlg::NoTranslation);
      auto stop = std::chrono::steady_clock::now();
      double timeWithoutTranslation = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolver.calclMultipoleExpansions(M2MAlg::ComplexTranslation);
      stop = std::chrono::steady_clock::now();
      double timeWithComplexTranslation = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolver.calclMultipoleExpansions(M2MAlg::RealTranslation);
      stop = std::chrono::steady_clock::now();
      double timeWithRealTranslation = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolver.calclMultipoleExpansions(M2MAlg::Layers, M2MDevice::CPU);
      stop = std::chrono::steady_clock::now();
      double timeWithLayersCPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolver.calclMultipoleExpansions(M2MAlg::Layers, M2MDevice::GPU);
      stop = std::chrono::steady_clock::now();
      double timeWithLayersGPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolver.calclMultipoleExpansions(M2MAlg::Matrices, M2MDevice::CPU);
      stop = std::chrono::steady_clock::now();
      double timeWithMatricesCPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolver.calclMultipoleExpansions(M2MAlg::Matrices, M2MDevice::GPU);
      stop = std::chrono::steady_clock::now();
      double timeWithMatricesGPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      //multipoleSolverWithMatricesAda.calclMultipoleExpansions(M2MAlg::Matrices, M2MDevice::Adaptive);
      stop = std::chrono::steady_clock::now();
      double timeWithMatricesAda = getTime(start, stop);

      std::cout << " " << octreeLeafCapacity;
      std::cout << " " << timeWithoutTranslation;
      std::cout << " " << timeWithComplexTranslation;
      std::cout << " " << timeWithRealTranslation;
      std::cout << " " << timeWithLayersCPU;
      std::cout << " " << timeWithLayersGPU;
      std::cout << " " << timeWithMatricesCPU;
      std::cout << " " << timeWithMatricesGPU << std::endl;
      //std::cout << " " << timeWithMatricesAda << std::endl;
   }
}


void calculationTimeForLocalMultipolesByNodeCount()
{
   constexpr double torusRadius = 2;
   constexpr double torusSectionWidth = 0.2;

   size_t w = 15;

   std::cout << std::setw(w) << "leaf capacity";
   std::cout << std::setw(w) << "w/t";
   std::cout << std::setw(w) << "Complex";
   std::cout << std::setw(w) << "real";
   std::cout << std::setw(w) << "layersCPU";
   std::cout << std::setw(w) << "layersGPU";
   std::cout << std::setw(w) << "matricesCPU";
   std::cout << std::setw(w) << "matricesGPU";
   //std::cout << std::setw(w) << "matricesAda";
   std::cout << std::endl;

   std::cout << std::fixed;

   for(size_t i = 1; i < 15; i++)
   {
      Torus torus(torusRadius, torusSectionWidth, 40 * i, 16, 16);
      auto bq = readBasisQuadratures();
      auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

      MultipoleSolver multipoleSolver(quadratures, 8);
      multipoleSolver.calcMultipoleExpansionsAtLeaves();
      multipoleSolver.log = false;

      auto start = std::chrono::steady_clock::now();
      //multipoleSolver.calclMultipoleExpansions(M2MAlg::NoTranslation);
      auto stop = std::chrono::steady_clock::now();
      double timeWithoutTranslation = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      //multipoleSolver.calclMultipoleExpansions(M2MAlg::ComplexTranslation);
      stop = std::chrono::steady_clock::now();
      double timeWithComplexTranslation = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      //multipoleSolver.calclMultipoleExpansions(M2MAlg::RealTranslation);
      stop = std::chrono::steady_clock::now();
      double timeWithRealTranslation = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolver.calclMultipoleExpansions(M2MAlg::Layers, M2MDevice::CPU);
      stop = std::chrono::steady_clock::now();
      double timeWithLayersCPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolver.calclMultipoleExpansions(M2MAlg::Layers, M2MDevice::GPU);
      stop = std::chrono::steady_clock::now();
      double timeWithLayersGPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolver.calclMultipoleExpansions(M2MAlg::Matrices, M2MDevice::CPU);
      stop = std::chrono::steady_clock::now();
      double timeWithMatricesCPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolver.calclMultipoleExpansions(M2MAlg::Matrices, M2MDevice::GPU);
      stop = std::chrono::steady_clock::now();
      double timeWithMatricesGPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      //multipoleSolverWithMatricesAda.calclMultipoleExpansions(M2MAlg::Matrices, M2MDevice::Adaptive);
      stop = std::chrono::steady_clock::now();
      double timeWithMatricesAda = getTime(start, stop);

      std::cout << " " << multipoleSolver.getOctreeNodeCount();
      std::cout << " " << timeWithoutTranslation;
      std::cout << " " << timeWithComplexTranslation;
      std::cout << " " << timeWithRealTranslation;
      std::cout << " " << timeWithLayersCPU;
      std::cout << " " << timeWithLayersGPU;
      std::cout << " " << timeWithMatricesCPU;
      std::cout << " " << timeWithMatricesGPU << std::endl;
      //std::cout << " " << timeWithMatricesAda << std::endl;
   }
}

std::pair<double, Vector3> wholeTimeForIntegrals(const std::vector<Vector3>& points, 
                           std::vector<Quadrature>& quadratures)
{
   Vector3 res;
   auto start = std::chrono::steady_clock::now();

   for(size_t p = 0; p < points.size(); p++)
   {
      res += math::calcBioSavartLaplace(current, points[p], quadratures);
   }

   auto stop = std::chrono::steady_clock::now();
   
   std::cout << res.x << " ";

   return std::pair<double, Vector3>
      (std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count() * 1e-9, res);
}

std::pair<double, Vector3>  timeForMultipoles(
   const std::vector<Vector3>& points,
   std::vector<Quadrature>& quadratures,
   M2MAlg alg, 
   M2MDevice device)
{
   Vector3 res;
   MultipoleSolver multipoleSolver(quadratures);
   multipoleSolver.calclMultipoleExpansions(alg, device);

   auto start = std::chrono::steady_clock::now();

   for(size_t p = 0; p < points.size(); p++)
   {
      res += multipoleSolver.calcB(current, points[p]);
   }

   auto stop = std::chrono::steady_clock::now();

   return { getTime(start, stop), res };
}

void NMResearch()
{
   const double torusRadius = 2;
   const double torusSectionWidth = 0.2;
   Vector3 begin(3, 1, 2);
   Vector3 end(0, 0, 0);
   BasisQuadratures bq = readBasisQuadratures();

   std::cout << " " << "NM";
   std::cout << " " << "w/t";
   std::cout << " " << "Complex";
   std::cout << " " << "real";
   std::cout << " " << "layersCPU";
   std::cout << " " << "layersGPU";
   std::cout << " " << "matricesCPU";
   std::cout << " " << "matricesGPU";
   std::cout << std::endl;

   for(size_t i = 5; i < 12; i++)
   {
      Torus torus(torusRadius, torusSectionWidth, pow(2, i + 1), 4, 4);
      auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

      int pointsCount = 4 * pow(2, i);
      auto points = createPoints(begin, end, pointsCount);

      std::cout << std::fixed;
      std::cout << std::setw(5) << pointsCount * quadratures.size() << " ";
      std::cout << std::scientific;
      std::cout << std::setw(8) << wholeTimeForIntegrals(points, quadratures).first << " ";
      std::cout << std::setw(8) << timeForMultipoles(points, quadratures, M2MAlg::NoTranslation, M2MDevice::CPU).first << " ";
      std::cout << std::setw(8) << timeForMultipoles(points, quadratures, M2MAlg::ComplexTranslation, M2MDevice::CPU).first << " ";
      std::cout << std::setw(8) << timeForMultipoles(points, quadratures, M2MAlg::RealTranslation, M2MDevice::CPU).first << " ";
      std::cout << std::setw(8) << timeForMultipoles(points, quadratures, M2MAlg::Layers, M2MDevice::CPU).first << " ";
      std::cout << std::setw(8) << timeForMultipoles(points, quadratures, M2MAlg::Layers, M2MDevice::GPU).first << " ";
      std::cout << std::setw(8) << timeForMultipoles(points, quadratures, M2MAlg::Matrices, M2MDevice::CPU).first << " ";
      std::cout << std::setw(8) << timeForMultipoles(points, quadratures, M2MAlg::Matrices, M2MDevice::GPU).first << std::endl;
   }
}

void layerCalculationsPrecision()
{
   Torus torus = createTorus();
   BasisQuadratures bq = readBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);
   MultipoleSolver multipoleSolverCPU(quadratures);
   MultipoleSolver multipoleSolverGPU(quadratures);

   Vector3 point(3, 1, 2);

   Vector3 byIntegration = math::calcBioSavartLaplace(current, point, quadratures);

   multipoleSolverCPU.calclMultipoleExpansions(M2MAlg::Layers, M2MDevice::CPU);
   Vector3 byMultipolesWithLayersCPU = multipoleSolverCPU.calcB(current, point);

   multipoleSolverGPU.calclMultipoleExpansions(M2MAlg::Layers, M2MDevice::GPU);
   Vector3 byMultipolesWithLayersGPU = multipoleSolverGPU.calcB(current, point);

   std::cout << std::setw(20) << "point ";
   point.printWithWidth(std::cout, 6);
   std::cout << std::scientific << std::endl;
   std::cout << std::setw(40) << "integration " << byIntegration << std::endl;
   std::cout << std::setw(40) << "multipoles with layers CPU" << byMultipolesWithLayersCPU << std::endl;
   std::cout << std::setw(40) << "multipoles with layers GPU" << byMultipolesWithLayersGPU << std::endl;
}

void layerCalculationTime()
{
   Torus torus = createTorus();
   BasisQuadratures bq = readBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

   Vector3 point(3, 1, 2);

   size_t w = 15;

   for(size_t i = 1; i < 2; i++)
   {
      int octreeLeafCapacity = pow(2, i);
      MultipoleSolver multipoleSolverCPU(quadratures, octreeLeafCapacity);
      MultipoleSolver multipoleSolverGPU(quadratures, octreeLeafCapacity);

      std::cout << std::setw(w) << "leaf capacity:";
      std::cout << std::setw(w) << octreeLeafCapacity << std::endl;

      auto start = std::chrono::steady_clock::now();
      multipoleSolverCPU.calclMultipoleExpansions(M2MAlg::Layers, M2MDevice::CPU);
      auto stop = std::chrono::steady_clock::now();
      double timeWithCPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolverGPU.calclMultipoleExpansions(M2MAlg::Layers, M2MDevice::GPU);
      stop = std::chrono::steady_clock::now();
      double timeWithGPU = getTime(start, stop);

      std::cout << std::setw(w) << "total time CPU:";
      std::cout << std::setw(w) << timeWithCPU << std::endl;
      std::cout << std::setw(w) << "total time GPU:";
      std::cout << std::setw(w) << timeWithGPU << std::endl;
   }
}

void matrixCalculationTime()
{
   Torus torus = createTorus();
   BasisQuadratures bq = readBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

   Vector3 point(3, 1, 2);

   size_t w = 15;

   for(size_t i = 1; i < 2; i++)
   {
      int octreeLeafCapacity = pow(2, i);
      MultipoleSolver multipoleSolverCPU(quadratures, octreeLeafCapacity);
      //MultipoleSolver multipoleSolverGPU(quadratures, quadratureOctreeLeafCapacity);
      //MultipoleSolver multipoleSolverAda(quadratures, quadratureOctreeLeafCapacity);

      std::cout << std::setw(w) << "leaf capacity:";
      std::cout << std::setw(w) << octreeLeafCapacity << std::endl;

      auto start = std::chrono::steady_clock::now();
      multipoleSolverCPU.calclMultipoleExpansions(M2MAlg::Matrices, M2MDevice::CPU);
      auto stop = std::chrono::steady_clock::now();
      double timeWithCPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      //multipoleSolverGPU.calclMultipoleExpansions(M2MAlg::Matrices, M2MDevice::GPU);
      stop = std::chrono::steady_clock::now();
      double timeWithGPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      //multipoleSolverAda.calclMultipoleExpansions(M2MAlg::Matrices, M2MDevice::Adaptive);
      stop = std::chrono::steady_clock::now();
      double timeWithAda = getTime(start, stop);

      std::cout << std::setw(w) << "total time CPU:";
      std::cout << std::setw(w) << timeWithCPU << std::endl;
      std::cout << std::setw(w) << "total time GPU:";
      std::cout << std::setw(w) << timeWithGPU << std::endl;
      std::cout << std::setw(w) << "total time Ada:";
      std::cout << std::setw(w) << timeWithAda << std::endl;
   }
}

void layerMatrixCalculationTime(M2MDevice device)
{
   Torus torus = createTorus();
   BasisQuadratures bq = readBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

   Vector3 point(3, 1, 2);

   size_t w = 15;

   for(size_t i = 2; i < 3; i++)
   {
      int octreeLeafCapacity = pow(2, i);
      MultipoleSolver multipoleSolverLayers(quadratures, octreeLeafCapacity);
      MultipoleSolver multipoleSolverMatrices(quadratures, octreeLeafCapacity);

      std::cout << std::setw(w) << "leaf capacity:";
      std::cout << std::setw(w) << octreeLeafCapacity << std::endl;

      auto start = std::chrono::steady_clock::now();
      multipoleSolverLayers.calclMultipoleExpansions(M2MAlg::Layers, device);
      auto stop = std::chrono::steady_clock::now();
      double timeWithLayers = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolverMatrices.calclMultipoleExpansions(M2MAlg::Matrices, device);
      stop = std::chrono::steady_clock::now();
      double timeWithMatrices = getTime(start, stop);

      
      std::cout << "device" << std::setw(w) << (device == M2MDevice::GPU ? "GPU" : "CPU") << std::endl;
      std::cout << "layers:  " << std::setw(w) << timeWithLayers << std::endl;
      std::cout << "matrices:" << std::setw(w) << timeWithMatrices << std::endl;
   }
}

void cudaAddingTest()
{
   std::vector<real> a = { 1, 2, 3, 4 };
   std::vector<real> b = { 10, 20, 30, 100 };
   std::cout << kernels::addVectors(a, b);
}

void printIndeces()
{
   int n = 10;
   int N = (n + 1) * (n + 1);
   std::vector<std::vector<bool>> matrix(N, std::vector<bool>(N));
   std::ofstream foutIdx("results/indeces.txt");
   std::ofstream foutMatrix("results/regular_matrix.txt");
   int w = 1;

   for(int l = 0; l <= n; l++)
   {
      for(int m = -l; m <= l; m++)
      {
         std::stringstream line;

         line << "M(" << std::setw(w) << l << "," << std::setw(w) << m << ") = ";

         for(int lambda = 0; lambda <= l; lambda++)
         {
            int dl = l - lambda;
            line << "[";

            for(int mu = -lambda; mu <= lambda; mu++)
            {
               int dm = m - mu;

               if(-dl <= dm && dm <= dl)
               {
                  line << "R(" << std::setw(w) << lambda << "," << std::setw(w) << mu << ")";
                  line << "M(" << std::setw(w) << dl << "," << std::setw(w) << dm << ")";
                  
                  if(mu != lambda && lambda != l)
                     line << " + ";
                  else
                     line << "] + ";

                  matrix[l * l + l + m][dl * dl + dl + dm] = true;
               }
            }
         }

         line << std::endl;
         std::string newStr = line.str();
         newStr.erase(newStr.end() - 3);
         line.str(newStr);
         foutIdx << line.str();
      }

      foutIdx << "------------------------------------" << std::endl;
   }

   //for(int l = 0; l <= _order; l++)
   //{
   //   for(int m = -l; m <= l; m++)
   //   {
   //      for(int lambda = 0; lambda <= l; lambda++)
   //      {
   //         for(int mu = -lambda; mu <= lambda; mu++)
   //         {
   //            foutMatrix << (matrix[l * l + l + m][lambda * lambda + lambda + mu] 
   //                           ? "O" : ".") << " ";
   //         }
   //      }
   //      foutMatrix << std::endl;
   //   }
   //}

   std::vector<int> indx;
   int id = 1;

   while(id < N)
   {
      indx.push_back((id * id) - 1);
      id++;
   }
   

   for(int i = 0; i < N; i++)
   {
      for(int j = 0; j < N; j++)
      {
         if(j <= i)
         {
            foutMatrix << (matrix[i][j] ? "O" : ".");

            for(size_t k = 0; k < indx.size(); k++)
            {
               if(indx[k] == j)
               {
                  foutMatrix << "|";
                  break;
               }
            }
         }
      }

      foutMatrix << std::endl;

      for(size_t k = 0; k < indx.size(); k++)
      {
         if(indx[k] == i)
         {
            for(size_t o = 0; o < N + n + 1; o++)
            {
               foutMatrix << "-";
            }

            foutMatrix << std::endl;
            break;
         }
      }
   }

}

void matrixCalculationsPrecision()
{
   Torus torus = createTorus();
   BasisQuadratures bq = readBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);
   MultipoleSolver multipoleSolver(quadratures);
   multipoleSolver.calcMultipoleExpansionsAtLeaves();
   //MultipoleSolver multipoleSolverAda(quadratures);

   Vector3 point(3, 1, 2);

   Vector3 byIntegration = math::calcBioSavartLaplace(current, point, quadratures);

   multipoleSolver.calclMultipoleExpansions(M2MAlg::Matrices, M2MDevice::CPU);
   Vector3 byMultipolesWithLayersCPU = multipoleSolver.calcB(current, point);

   multipoleSolver.calclMultipoleExpansions(M2MAlg::Matrices, M2MDevice::GPU);
   Vector3 byMultipolesWithLayersGPU = multipoleSolver.calcB(current, point);

   //multipoleSolverAda.calclMultipoleExpansions(M2MAlg::Matrices, M2MDevice::Adaptive);
   //Vector3 byMultipolesWithLayersAda = multipoleSolverAda.calcB(current, point);

   std::cout << std::setw(20) << "point ";
   point.printWithWidth(std::cout, 6);
   std::cout << std::scientific << std::endl;
   std::cout << std::setw(40) << "integration " << byIntegration << std::endl;
   std::cout << std::setw(40) << "multipoles with matrices CPU" << byMultipolesWithLayersCPU << std::endl;
   std::cout << std::setw(40) << "multipoles with matrices GPU" << byMultipolesWithLayersGPU << std::endl;
   //std::cout << std::setw(40) << "multipoles with matrices Ada" << byMultipolesWithLayersAda << std::endl;
}

void translationTest()
{
   Vector3 a(3, 2, 1);
   Vector3 b(1, 7, 4);
   
   auto seriesA = Harmonics(10, a).sphericalHarmonics();
   auto complexSeriesA1 = Harmonics::realToComplex(seriesA);
   auto realSeriesA1 = Harmonics::complexToReal(complexSeriesA1);
   
   std::vector<Complex> temp1(121);
   cblas_scopy(121, seriesA.data().data(), 1, (float*)temp1.data(), 2);
   auto realToComplexMatrix2D = Harmonics::calcRealToComplexMatrix2D(10);
   auto complexSeriesA2 = math::mult(temp1, realToComplexMatrix2D);
   auto complexToRealMatrix2D = Harmonics::calcComplexToRealMatrix2D(10);
   auto temp2 = math::mult(complexSeriesA2, complexToRealMatrix2D);
   std::vector<real> realSeriesA2(121);
   cblas_scopy(121, (float*)temp2.data(), 2, realSeriesA2.data(), 1);

   auto realToComplexMatrixTransposed1D =
      Harmonics::calcRealToComplexTransitionMatrix1D(10);
   auto complexToRealMatrixTransposed1D =
      Harmonics::calcComplexToRealTransitionMatrix1D(10);

   auto matricesMultiplied = math::multMatricesAsVectors(
      complexToRealMatrixTransposed1D,
      realToComplexMatrixTransposed1D,
      121, 121, 121
   );

   auto deb = math::multMatricesAsVectors(
      temp1,
      matricesMultiplied,
      121, 1, 121);
}

void multipoleToLocalTest()
{
   int n = 10;
   Vector3 point1(3, 1, 2);
   Vector3 point2(4, 5, 1);

   auto r1 = Harmonics::calcRegularSolidHarmonics(n, point1);
   auto c1 = Harmonics::realToComplex(r1);

   auto r2 = Harmonics::calcRegularSolidHarmonics(n, point2);
   auto c2 = Harmonics::realToComplex(r2);

   auto t = Harmonics::complexToReal(
      MultipoleTranslator::translateLocal(c1, c2));
}

void FMMPrecisionTest()
{
   Torus torus = createTorus();
   BasisQuadratures bq = readBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);
   //Vector3 begin(10, 10, 10);
   //Vector3 end(9, 9, 9);
   Vector3 begin(4, 4, 4);
   Vector3 end(1, 1, 1);
   //auto points = createPoints(begin, end, 10);
   //std::vector<Vector3> points = {{2, 2, 0}};
   //std::vector<Vector3> points = {{3, 3, 3}};
   auto points = createRandomPoints(Box({ 0, 0, 0 }, { 2, 2, 2 }), 128);
   //auto points = createRandomPoints(Box({ 10, 10, 10 }, { 2, 2, 2 }), 100);
   //auto points = std::vector<Vector3>({ Vector3(10.5, 5, 8) });

   FastMultipoleSolver multipoleSolver(quadratures, points, 64, 32);
   multipoleSolver.log = false;
   multipoleSolver.calcMultipoleExpansionsAtLeaves();
   multipoleSolver.calclMultipoleExpansions(M2MAlg::Matrices, M2MDevice::GPU);
   multipoleSolver.calclLocalMultipoleExpansions(M2LAlg::ComplexTranslation, M2MDevice::CPU);

   auto fmmResults = multipoleSolver.calcA(current);
   
   real averageAbsoluteError = 0;
   real averageRelativeError = 0;

   for(size_t i = 0; i < points.size(); ++i)
   {
      auto point = fmmResults[i].first;
      auto aInPointByFmm = fmmResults[i].second;
      auto byMultipolesWithMatricesGPU = multipoleSolver.calcA(current, point);

      real absoluteError = (aInPointByFmm - byMultipolesWithMatricesGPU).length();
      real relativeError = 100 * absoluteError / byMultipolesWithMatricesGPU.length();

      averageAbsoluteError += absoluteError;
      averageRelativeError += relativeError;

      test::printSeparateLine(std::cout, 100);
      std::cout << std::scientific;
      std::cout << std::setw(40) << "point: " << point << "\n";
      std::cout << std::setw(40) << "multipoles with matrices GPU: " << byMultipolesWithMatricesGPU << "\n";
      std::cout << std::setw(40) << "fmm: " << aInPointByFmm << " ";
      std::cout << std::setw(10) << absoluteError << "\n";
   }

   std::cout << "averageAbsoluteError: " << averageAbsoluteError / points.size() << "\n";
   std::cout << "averageRelativeError: " << averageRelativeError / points.size() << "\n";
}

void FFMTimeTest()
{
   Torus torus = createTorus();
   BasisQuadratures bq = readBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

   std::cout << std::setw(5);
   std::cout << std::setw(16) << "difference";
   std::cout << std::setw(16) << "noFMMTime" << std::setw(16) << "FMMTime" << "\n";
   test::printSeparateLine(std::cout, 70);

   for(size_t i = 0; i < 33; i++)
   {
      size_t pointCount = pow(2, i);
      //size_t pointCount = (i + 1) * quadratures.size();
      auto points = createRandomPoints(Box({ 0, 0, 0 }, { 2, 2, 2 }), pointCount);
      FastMultipoleSolver multipoleSolver(quadratures, points, 128, 32);
      multipoleSolver.log = false;

      auto commonPartStart = std::chrono::steady_clock::now();
      multipoleSolver.calcMultipoleExpansionsAtLeaves();
      multipoleSolver.calclMultipoleExpansions(M2MAlg::Matrices, M2MDevice::GPU);
      auto commonPartStop = std::chrono::steady_clock::now();
      double commonPartTime = test::getTime(commonPartStart, commonPartStop);

      auto fmmPartStart = std::chrono::steady_clock::now();
      multipoleSolver.calclLocalMultipoleExpansions(M2LAlg::ComplexTranslation, M2MDevice::CPU);
      auto fmmPartStop = std::chrono::steady_clock::now();
      double fmmPartTime = test::getTime(fmmPartStart, fmmPartStop);


      Vector3 noFMMRes;
      Vector3 FMMRes;

      auto start = std::chrono::steady_clock::now();
      /*for(auto& point : points)
      {
         noFMMRes += multipoleSolver.calcA(current, point);
      }*/
      auto stop = std::chrono::steady_clock::now();
      double noFMMTime = test::getTime(start, stop);

      start = std::chrono::steady_clock::now();
      auto fmmResults = multipoleSolver.calcA(current);
      stop = std::chrono::steady_clock::now();
      double FMMTime = test::getTime(start, stop);

      for (auto &[point, fmmResult] : fmmResults)
      {
         FMMRes += fmmResult;
      }

      std::cout << std::fixed << std::setw(9) << pointCount;
      std::cout << std::scientific;
      //std::cout << std::setw(16) << 100 * (noFMMRes - FMMRes).length() / noFMMRes.length();
      //std::cout << std::setw(16) << noFMMTime;
      //std::cout << std::setw(16) << FMMTime + fmmPartTime << "\n";
      std::cout << "\t" << FMMTime + fmmPartTime << "\n";
   }
}

int main()
{
   //NMResearch();
   //timeResearchForMorePoints();
   //comparisonToTelmaIntegrals();
   //octreeFormingTime();
   //calculationTimeForMultipolesInLeaves();
   //comparisonBetweenMethodsOnPrecision();
   //calculationTimeForLocalMultipolesByNodeCount();
   //layerCalculationsPrecision();
   //matrixCalculationsPrecision();

   //multipoleToLocalTest();

   //layerCalculationTime();
   //matrixCalculationTime();
   //layerMatrixCalculationTime(M2MDevice::CPU);
   //layerMatrixCalculationTime(M2MDevice::GPU);
   //compareWithMatrixMultiplication();

   FMMPrecisionTest();
   //FFMTimeTest();
}
