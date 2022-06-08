#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>

#include "cblas.h"
#include "cylinder.hpp"
#include "exeptions.hpp"
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
   BasisQuadratures bq = readTetrahedronBasisQuadratures();
   auto telmaResults = readTelmaResults("results/telma_results.txt");
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

   std::vector<Vector3> initialPoints(telmaResults.size());

   for(size_t i = 0; i < telmaResults.size(); i++)
   {
      initialPoints[i] = telmaResults[i].first;
   }

   MultipoleSolver multipoleSolverNoExp(quadratures);
   MultipoleSolver multipoleSolverComplex(quadratures);
   MultipoleSolver multipoleSolverReal(quadratures);
   MultipoleSolver multipoleSolverLayersCPU(quadratures);
   MultipoleSolver multipoleSolverLayersGPU(quadratures);
   MultipoleSolver multipoleSolverMatricesCPU(quadratures);
   MultipoleSolver multipoleSolverMatricesGPU(quadratures);

   FastMultipoleSolver fmmSolver(quadratures, initialPoints, Problem::BioSavartLaplace, 1000, 32);

   multipoleSolverNoExp.calclMultipoleExpansions(M2MAlg::NoTranslation);
   multipoleSolverComplex.calclMultipoleExpansions(M2MAlg::ComplexTranslation);
   multipoleSolverReal.calclMultipoleExpansions(M2MAlg::RealTranslation);
   multipoleSolverLayersCPU.calclMultipoleExpansions(M2MAlg::Layers, Device::CPU);
   multipoleSolverLayersGPU.calclMultipoleExpansions(M2MAlg::Layers, Device::GPU);
   multipoleSolverMatricesCPU.calclMultipoleExpansions(M2MAlg::Matrices, Device::CPU);
   multipoleSolverMatricesGPU.calclMultipoleExpansions(M2MAlg::Matrices, Device::GPU);

   fmmSolver.calclMultipoleExpansions(M2MAlg::Matrices, Device::GPU);
   fmmSolver.calcLocalMultipoleExpansions(M2LAlg::ComplexTranslation, Device::CPU);

   auto points = fmmSolver.calcB(current);

   real sumErrorIntegration = 0;
   real sumErrorNoTranslation = 0;
   real sumErrorComplexTranslation = 0;
   real sumErrorRealTranslation = 0;
   real sumErrorLayersTranslationCPU = 0;
   real sumErrorLayersTranslationGPU = 0;
   real sumErrorMatricesTranslationCPU = 0;
   real sumErrorMatricesTranslationGPU = 0;
   real sumFmmSolver = 0;

   size_t n = telmaResults.size();

   for(size_t i = 0; i < n; i++)
   {
      auto point = points[i].first;
      Vector3 telmaB;

      for(size_t j = 0; j < n; j++)
      {
         if(point.x == telmaResults[j].first.x &&
            point.y == telmaResults[j].first.y &&
            point.z == telmaResults[j].first.z)
         {
            telmaB = telmaResults[j].second * math::MU0;
            break;
         }
      }

      /*auto point = telmaResults[i].first;
      auto telmaB = telmaResults[i].second * math::MU0;*/

      Vector3 integration = math::calcBioSavartLaplace(current, point, quadratures);
      Vector3 noTranslation = multipoleSolverNoExp.calcB(current, point);
      Vector3 complexTranslation = multipoleSolverComplex.calcB(current, point);
      Vector3 realTranslation = multipoleSolverReal.calcB(current, point);
      Vector3 layersTranslationCPU = multipoleSolverLayersCPU.calcB(current, point);
      Vector3 layersTranslationGPU = multipoleSolverLayersGPU.calcB(current, point);
      Vector3 matricesTranslationCPU = multipoleSolverMatricesCPU.calcB(current, point);
      Vector3 matricesTranslationGPU = multipoleSolverMatricesGPU.calcB(current, point);

      real errorIntegration = 100 * (integration - telmaB).length() / telmaB.length();
      real errorNoTranslation = 100 * (noTranslation - telmaB).length() / telmaB.length();
      real errorComplexTranslation = 100 * (complexTranslation - telmaB).length() / telmaB.length();
      real errorRealTranslation = 100 * (realTranslation - telmaB).length() / telmaB.length();
      real errorLayersTranslationCPU = 100 * (layersTranslationCPU - telmaB).length() / telmaB.length();
      real errorLayersTranslationGPU = 100 * (layersTranslationGPU - telmaB).length() / telmaB.length();
      real errorMatricesTranslationCPU = 100 * (matricesTranslationCPU - telmaB).length() / telmaB.length();
      real errorMatricesTranslationGPU = 100 * (matricesTranslationGPU - telmaB).length() / telmaB.length();
      real errorFmm = 100 * (points[i].second - telmaB).length() / telmaB.length();

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
      std::cout << std::setw(16) << errorFmm;
      std::cout << std::endl;

      sumErrorIntegration += errorIntegration;
      sumErrorNoTranslation += errorNoTranslation;
      sumErrorComplexTranslation += errorComplexTranslation;
      sumErrorRealTranslation += errorRealTranslation;
      sumErrorLayersTranslationCPU += errorLayersTranslationCPU;
      sumErrorLayersTranslationGPU += errorLayersTranslationGPU;
      sumErrorMatricesTranslationCPU += errorMatricesTranslationCPU;
      sumErrorMatricesTranslationGPU += errorMatricesTranslationGPU;
      sumFmmSolver += errorFmm;
   }

   std::cout << sumErrorIntegration / n << std::endl;
   std::cout << sumErrorNoTranslation / n << std::endl;
   std::cout << sumErrorComplexTranslation / n << std::endl;
   std::cout << sumErrorRealTranslation / n << std::endl;
   std::cout << sumErrorLayersTranslationCPU / n << std::endl;
   std::cout << sumErrorLayersTranslationGPU / n << std::endl;
   std::cout << sumErrorMatricesTranslationCPU / n << std::endl;
   std::cout << sumErrorMatricesTranslationGPU / n << std::endl;
   std::cout << sumFmmSolver / n << std::endl;
}

void comparisonBetweenMethodsOnPrecision()
{
   Torus torus = createTorus();
   BasisQuadratures bq = readTetrahedronBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);
   MultipoleSolver multipoleSolver(quadratures, Problem::BioSavartLaplace, 10);
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

   multipoleSolver.calclMultipoleExpansions(M2MAlg::Layers, Device::GPU);
   Vector3 byMultipolesWithLayersGPU = multipoleSolver.calcB(current, point);

   multipoleSolver.calclMultipoleExpansions(M2MAlg::Matrices, Device::CPU);
   Vector3 byMultipolesWithMatricesCPU = multipoleSolver.calcB(current, point);

   multipoleSolver.calclMultipoleExpansions(M2MAlg::Matrices, Device::GPU);
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
      auto bq = readTetrahedronBasisQuadratures();
      auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

      auto start = std::chrono::steady_clock::now();
      MultipoleSolver multipoleSolver(quadratures, Problem::BioSavartLaplace, 1000);
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
      auto bq = readTetrahedronBasisQuadratures();
      auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

      MultipoleSolver multipoleSolver(quadratures, Problem::BioSavartLaplace, 16);

      auto start = std::chrono::steady_clock::now();
      multipoleSolver.calcMultipoleExpansionsAtLeaves();
      auto stop = std::chrono::steady_clock::now();
      double timeForMultipolesInLeaves = getTime(start, stop);

      std::cout << multipoleSolver.getQuadratureOctreeNodeCount() << " ";
      std::cout << timeForMultipolesInLeaves << std::endl;
   }
}

void calculationTimeForLocalMultipolesByLeafCapacity()
{
   auto torus = createTorus();
   auto bq = readTetrahedronBasisQuadratures();
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
      MultipoleSolver multipoleSolver(quadratures, Problem::BioSavartLaplace, octreeLeafCapacity);
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
      multipoleSolver.calclMultipoleExpansions(M2MAlg::Layers, Device::CPU);
      stop = std::chrono::steady_clock::now();
      double timeWithLayersCPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolver.calclMultipoleExpansions(M2MAlg::Layers, Device::GPU);
      stop = std::chrono::steady_clock::now();
      double timeWithLayersGPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolver.calclMultipoleExpansions(M2MAlg::Matrices, Device::CPU);
      stop = std::chrono::steady_clock::now();
      double timeWithMatricesCPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolver.calclMultipoleExpansions(M2MAlg::Matrices, Device::GPU);
      stop = std::chrono::steady_clock::now();
      double timeWithMatricesGPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      //multipoleSolverWithMatricesAda.calclMultipoleExpansions(M2MAlg::Matrices, Device::Adaptive);
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
      auto bq = readTetrahedronBasisQuadratures();
      auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

      MultipoleSolver multipoleSolver(quadratures, Problem::BioSavartLaplace, 8);
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
      multipoleSolver.calclMultipoleExpansions(M2MAlg::Layers, Device::CPU);
      stop = std::chrono::steady_clock::now();
      double timeWithLayersCPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolver.calclMultipoleExpansions(M2MAlg::Layers, Device::GPU);
      stop = std::chrono::steady_clock::now();
      double timeWithLayersGPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolver.calclMultipoleExpansions(M2MAlg::Matrices, Device::CPU);
      stop = std::chrono::steady_clock::now();
      double timeWithMatricesCPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolver.calclMultipoleExpansions(M2MAlg::Matrices, Device::GPU);
      stop = std::chrono::steady_clock::now();
      double timeWithMatricesGPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      //multipoleSolverWithMatricesAda.calclMultipoleExpansions(M2MAlg::Matrices, Device::Adaptive);
      stop = std::chrono::steady_clock::now();
      double timeWithMatricesAda = getTime(start, stop);

      std::cout << " " << multipoleSolver.getQuadratureOctreeNodeCount();
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
   
   return std::pair<double, Vector3>
      (std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count() * 1e-9, res);
}

std::pair<double, Vector3>  timeForMultipoles(
   const std::vector<Vector3>& points,
   std::vector<Quadrature>& quadratures,
   M2MAlg alg, 
   Device device)
{
   Vector3 res;
   auto start = std::chrono::steady_clock::now();
   MultipoleSolver multipoleSolver(quadratures);
   multipoleSolver.calclMultipoleExpansions(alg, device);

   for(size_t p = 0; p < points.size(); p++)
   {
      res += multipoleSolver.calcB(current, points[p]);
   }

   auto stop = std::chrono::steady_clock::now();

   return { getTime(start, stop), res };
}

void NMResearch1()
{
   const double torusRadius = 2;
   const double torusSectionWidth = 0.2;
   Vector3 begin(3, 1, 2);
   Vector3 end(0, 0, 0);
   BasisQuadratures bq = readTetrahedronBasisQuadratures();

   std::cout << std::setw(16) << "NM";
   std::cout << std::setw(16) << "w/t";
   std::cout << std::setw(16) << "matrices";
   std::cout << std::endl;

   double prevTimeIntegrals = 0;
   double prevTimeMatrices = 0;

   for(size_t i = 1; i < 15; i++)
   {
      Torus torus(torusRadius, torusSectionWidth, pow(2, i), 4, 4);
      auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

      int pointsCount = 3840 * pow(2, i - 1);
      auto points = createPoints(begin, end, pointsCount);

      std::cout << std::fixed;
      std::cout << std::setw(16) << pointsCount * quadratures.size() << " ";
      std::cout << std::scientific;

      if(i < 5)
         prevTimeIntegrals = wholeTimeForIntegrals(points, quadratures).first;
      else
         prevTimeIntegrals *= 4;

      if(i < 8)
         prevTimeMatrices = timeForMultipoles(points, quadratures, M2MAlg::Matrices, Device::GPU).first;
      else
         prevTimeMatrices = 0;

      std::cout << std::setw(16) << prevTimeIntegrals;
      std::cout << std::setw(16) << prevTimeMatrices << std::endl;
   }
}

void NMResearch2()
{
   const double torusRadius = 2;
   const double torusSectionWidth = 0.2;
   BasisQuadratures bq = readTetrahedronBasisQuadratures();

   std::cout << std::setw(16) << "NM";
   std::cout << std::setw(16) << "w/t";
   std::cout << std::setw(16) << "matrices";
   std::cout << std::setw(16) << "Fmm";
   std::cout << std::endl;

   double prevTimeIntegrals = 0;
   double prevTimeMatrices = 0;

   for(size_t i = 1; i < 15; i++)
   {
      Torus torus(torusRadius, torusSectionWidth, pow(2, i), 4, 4);
      auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

      int pointCount = 3840 * pow(2, i - 1);
      auto points = createRandomPoints(Box({ 0, 0, 0 }, { 2, 2, 2 }), pointCount);

      std::cout << std::fixed;
      std::cout << std::setw(16) << pointCount * quadratures.size() << " ";
      std::cout << std::scientific;

      /*if(i < 5)
         prevTimeIntegrals = wholeTimeForIntegrals(points, quadratures).first;
      else
         prevTimeIntegrals *= 4;*/
      

      if(i < 15)
         prevTimeMatrices = timeForMultipoles(points, quadratures, M2MAlg::ComplexTranslation, Device::GPU).first;
      else
         prevTimeMatrices *= 2;

      std::cout << std::setw(16) << prevTimeIntegrals;
      std::cout << std::setw(16) << prevTimeMatrices;

      auto start = std::chrono::steady_clock::now();
      FastMultipoleSolver fmmSolver(quadratures, points, Problem::BioSavartLaplace, 1000, 100);
      fmmSolver.calcMultipoleExpansionsAtLeaves();
      fmmSolver.calclMultipoleExpansions(M2MAlg::ComplexTranslation, Device::GPU);
      fmmSolver.calcLocalMultipoleExpansions(M2LAlg::ComplexTranslation, Device::CPU);
      auto fmmResults = fmmSolver.calcB(current);
      auto stop = std::chrono::steady_clock::now();
      auto time = getTime(start, stop);

      std::cout << std::setw(16) << time << std::endl;
   }
}

void layerCalculationsPrecision()
{
   Torus torus = createTorus();
   BasisQuadratures bq = readTetrahedronBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);
   MultipoleSolver multipoleSolverCPU(quadratures);
   MultipoleSolver multipoleSolverGPU(quadratures);

   Vector3 point(3, 1, 2);

   Vector3 byIntegration = math::calcBioSavartLaplace(current, point, quadratures);

   multipoleSolverCPU.calclMultipoleExpansions(M2MAlg::Layers, Device::CPU);
   Vector3 byMultipolesWithLayersCPU = multipoleSolverCPU.calcB(current, point);

   multipoleSolverGPU.calclMultipoleExpansions(M2MAlg::Layers, Device::GPU);
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
   BasisQuadratures bq = readTetrahedronBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

   Vector3 point(3, 1, 2);

   size_t w = 15;

   for(size_t i = 1; i < 2; i++)
   {
      int octreeLeafCapacity = pow(2, i);
      MultipoleSolver multipoleSolverCPU(quadratures, Problem::BioSavartLaplace, octreeLeafCapacity);
      MultipoleSolver multipoleSolverGPU(quadratures, Problem::BioSavartLaplace, octreeLeafCapacity);

      std::cout << std::setw(w) << "leaf capacity:";
      std::cout << std::setw(w) << octreeLeafCapacity << std::endl;

      auto start = std::chrono::steady_clock::now();
      multipoleSolverCPU.calclMultipoleExpansions(M2MAlg::Layers, Device::CPU);
      auto stop = std::chrono::steady_clock::now();
      double timeWithCPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolverGPU.calclMultipoleExpansions(M2MAlg::Layers, Device::GPU);
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
   BasisQuadratures bq = readTetrahedronBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

   Vector3 point(3, 1, 2);

   size_t w = 15;

   for(size_t i = 1; i < 2; i++)
   {
      int octreeLeafCapacity = pow(2, i);
      MultipoleSolver multipoleSolverCPU(quadratures, Problem::BioSavartLaplace, octreeLeafCapacity);
      //MultipoleSolver multipoleSolverGPU(quadratures, quadratureOctreeLeafCapacity);
      //MultipoleSolver multipoleSolverAda(quadratures, quadratureOctreeLeafCapacity);

      std::cout << std::setw(w) << "leaf capacity:";
      std::cout << std::setw(w) << octreeLeafCapacity << std::endl;

      auto start = std::chrono::steady_clock::now();
      multipoleSolverCPU.calclMultipoleExpansions(M2MAlg::Matrices, Device::CPU);
      auto stop = std::chrono::steady_clock::now();
      double timeWithCPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      //multipoleSolverGPU.calclMultipoleExpansions(M2MAlg::Matrices, Device::GPU);
      stop = std::chrono::steady_clock::now();
      double timeWithGPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      //multipoleSolverAda.calclMultipoleExpansions(M2MAlg::Matrices, Device::Adaptive);
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

void layerMatrixCalculationTime(Device device)
{
   Torus torus = createTorus();
   BasisQuadratures bq = readTetrahedronBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

   Vector3 point(3, 1, 2);

   size_t w = 15;

   for(size_t i = 2; i < 3; i++)
   {
      int octreeLeafCapacity = pow(2, i);
      MultipoleSolver multipoleSolverLayers(quadratures, Problem::BioSavartLaplace, octreeLeafCapacity);
      MultipoleSolver multipoleSolverMatrices(quadratures, Problem::BioSavartLaplace, octreeLeafCapacity);

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

      
      std::cout << "device" << std::setw(w) << (device == Device::GPU ? "GPU" : "CPU") << std::endl;
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
   BasisQuadratures bq = readTetrahedronBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);
   MultipoleSolver multipoleSolver(quadratures);
   multipoleSolver.calcMultipoleExpansionsAtLeaves();
   //MultipoleSolver multipoleSolverAda(quadratures);

   Vector3 point(3, 1, 2);

   Vector3 byIntegration = math::calcBioSavartLaplace(current, point, quadratures);

   multipoleSolver.calclMultipoleExpansions(M2MAlg::Matrices, Device::CPU);
   Vector3 byMultipolesWithLayersCPU = multipoleSolver.calcB(current, point);

   multipoleSolver.calclMultipoleExpansions(M2MAlg::Matrices, Device::GPU);
   Vector3 byMultipolesWithLayersGPU = multipoleSolver.calcB(current, point);

   //multipoleSolverAda.calclMultipoleExpansions(M2MAlg::Matrices, Device::Adaptive);
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
   /*Vector3 a(3, 2, 1);
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
      121, 1, 121);*/

   Vector3 a(3, 2, 1);
   Vector3 b(1, 7, 4);

   Vector3 translation = a - b;

   auto seriesA = Harmonics(10, a).sphericalHarmonics();
   auto reg = Harmonics::calcRegularSolidHarmonics(10, translation);
   
   auto seriesA1Translated = MultipoleTranslator::translateMultipole(
      seriesA, reg);

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
   BasisQuadratures bq = readTetrahedronBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);
   //Vector3 begin(10, 10, 10);
   //Vector3 end(9, 9, 9);
   Vector3 begin(4, 4, 4);
   Vector3 end(1, 1, 1);
   //auto points = createPoints(begin, end, 10);
   //std::vector<Vector3> points = {{2, 2, 0}};
   //std::vector<Vector3> points = {{3, 3, 3}};
   //auto points = createRandomPoints(Box({ 0, 0, 0 }, { 2, 2, 2 }), 256);
   auto points = createRandomPoints(Box({ 10, 10, 10 }, { 2, 2, 2 }), 256);
   //auto points = std::vector<Vector3>({ Vector3(10.5, 5, 8) });

   FastMultipoleSolver multipoleSolver(quadratures, points, Problem::BioSavartLaplace, 64, 32);
   multipoleSolver.log = false;
   multipoleSolver.calcMultipoleExpansionsAtLeaves();
   multipoleSolver.calclMultipoleExpansions(M2MAlg::Matrices, Device::GPU);
   //multipoleSolver.calcLocalMultipoleExpansions(M2LAlg::ComplexTranslation, Device::CPU);
   multipoleSolver.calcLocalMultipoleExpansions(M2LAlg::Matrices, Device::CPU);

   auto fmmResults = multipoleSolver.calcB(current);
   
   real averageAbsoluteError = 0;
   real averageRelativeError = 0;

   for(size_t i = 0; i < points.size(); ++i)
   {
      auto point = fmmResults[i].first;
      auto aInPointByFmm = fmmResults[i].second;
      auto byMultipolesWithMatricesGPU = multipoleSolver.calcB(current, point);

      real absoluteError = (aInPointByFmm - byMultipolesWithMatricesGPU).length();
      real relativeError = 100 * absoluteError / byMultipolesWithMatricesGPU.length();

      averageAbsoluteError += absoluteError;
      averageRelativeError += relativeError;

      test::printSeparateLine(std::cout, 100);
      std::cout << std::scientific;
      std::cout << std::setw(40) << "point: " << point << "\n";
      std::cout << std::setw(40) << "multipoles with matrices GPU: " << byMultipolesWithMatricesGPU << "\n";
      std::cout << std::setw(40) << "fmm: " << aInPointByFmm << " ";
      std::cout << std::setw(10) << relativeError << "\n";
   }

   std::cout << "averageAbsoluteError: " << averageAbsoluteError / points.size() << "\n";
   std::cout << "averageRelativeError: " << averageRelativeError / points.size() << "\n";
}

void FFMTimeTest()
{
   Torus torus = createTorus();
   BasisQuadratures bq = readTetrahedronBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

   std::cout << std::setw(9) << "points";
   std::cout << std::setw(16) << "errorComplex";
   //std::cout << std::setw(16) << "errorMatrices";
   std::cout << std::setw(16) << "noFMMTime";
   std::cout << std::setw(16) << "timeComplex" << "\n";
   //std::cout << std::setw(16) << "timeMatrices" << "\n";
   test::printSeparateLine(std::cout, 100);

   for(size_t i = 0; i < 33; i++)
   {
      size_t pointCount = pow(2, i);
      //size_t pointCount = (i + 1) * quadratures.size();
      //auto points = createRandomPoints(Box({ 0, 0, 0 }, { 2, 2, 2 }), pointCount);
      auto points = createRandomPoints(Box({ 0, 0, 0 }, { 2, 2, 2 }), pointCount);
      FastMultipoleSolver multipoleSolverComplex(quadratures, points, Problem::BioSavartLaplace, 128, 32);
      //FastMultipoleSolver multipoleSolverMatricesGPU(quadratures, points, 128, 32);
      multipoleSolverComplex.log = false;
      //multipoleSolverMatricesGPU.log = false;
      
      multipoleSolverComplex.calcMultipoleExpansionsAtLeaves();
      multipoleSolverComplex.calclMultipoleExpansions(M2MAlg::Matrices, Device::GPU);

      auto start = std::chrono::steady_clock::now();
      multipoleSolverComplex.calcLocalMultipoleExpansions(M2LAlg::ComplexTranslation, Device::GPU);
      auto stop = std::chrono::steady_clock::now();
      double fmmPartTimeComplex = test::getTime(start, stop);
      
      //multipoleSolverMatricesGPU.calcMultipoleExpansionsAtLeaves();
      //multipoleSolverMatricesGPU.calclMultipoleExpansions(M2MAlg::Matrices, Device::GPU);

      //start = std::chrono::steady_clock::now();
      //multipoleSolverMatricesGPU.calcLocalMultipoleExpansions(M2LAlg::Matrices, Device::GPU);
      //stop = std::chrono::steady_clock::now();
      //double fmmPartTimeMatrices = test::getTime(start, stop);

      Vector3 noFMMResComplex;
      Vector3 FMMResComplex;
      //Vector3 FMMResMatrices;

      start = std::chrono::steady_clock::now();
      for(auto& point : points)
      {
         noFMMResComplex += multipoleSolverComplex.calcB(current, point);
      }
      stop = std::chrono::steady_clock::now();
      double noFMMTime = test::getTime(start, stop);

      start = std::chrono::steady_clock::now();
      auto fmmComplexResults = multipoleSolverComplex.calcB(current);
      stop = std::chrono::steady_clock::now();
      double FMMTime = test::getTime(start, stop);

      //start = std::chrono::steady_clock::now();
      //auto fmmMatricesResults = multipoleSolverMatricesGPU.calcB(current);
      //stop = std::chrono::steady_clock::now();
      //double FMMTimeMatrices = test::getTime(start, stop);

      for (auto &[point, fmmResult] : fmmComplexResults)
      {
         FMMResComplex += fmmResult;
      }

      //for(auto& [point, fmmResult] : fmmMatricesResults)
      //{
      //   FMMResMatrices += fmmResult;
      //}

      std::cout << std::fixed << std::setw(9) << pointCount;
      std::cout << std::scientific;
      std::cout << std::setw(16) << 100 * (noFMMResComplex - FMMResComplex).length() / noFMMResComplex.length();
      //std::cout << std::setw(16) << 100 * (noFMMResComplex - FMMResMatrices).length() / noFMMResComplex.length();
      std::cout << std::setw(16) << noFMMTime;
      std::cout << std::setw(16) << FMMTime + fmmPartTimeComplex << "\n";
      //std::cout << std::setw(16) << FMMTimeMatrices + fmmPartTimeMatrices << "\n";
   }
}

void translationTest2()
{
   Vector3 a(4.2, 6.9, 2.1);
   Vector3 b(4, 7, 2);
   size_t order = 10;
   size_t harmonicLength = 121;
   size_t matrixElemCount = harmonicLength * harmonicLength;

   //auto trans = a - b;
   auto trans = a - b;
   auto harm = Harmonics::calcRegularSolidHarmonics(order, b);
   auto harmAComplex = Harmonics::realToComplex(harm);
   auto regular = Harmonics::realToComplex(Harmonics::calcRegularSolidHarmonics(order, trans));
   //auto realRegular = Harmonics::calcRegularSolidHarmonics(order, trans);
   //auto trueResult = MultipoleTranslator::translateMultipole(realRegular, harm).data();

   auto trueResult1 = Harmonics::complexToReal(
      MultipoleTranslator::translateLocal(harmAComplex, regular)).data();

   std::vector<Complex> regularMatrix(harmonicLength * harmonicLength);

   /*for(int l = 0; l <= order; l++)
   {
      for(int m = -l; m <= l; m++)
      {
         for(int lambda = 0; lambda <= l; lambda++)
         {
            int dl = l - lambda;

            for(int mu = -lambda; mu <= lambda; mu++)
            {
               int dm = m - mu;

               if(-dl <= dm && dm <= dl)
               {
                  regularMatrix[(l * l + l + m) + (dl * dl + dl + dm) * harmonicLength] =
                     regular.getHarmonic(lambda * lambda + lambda + mu) *
                     MultipoleTranslator::multipoleTranslationFactor(m, mu);
               }
            }
         }
      }
   }*/

   for(int l = 0; l <= regular.order(); l++)
   {
      for(int m = -l; m <= l; m++)
      {
         for(int lambda = 0; lambda <= l; lambda++)
         {
            int dl = lambda - l;
            if(dl >= 0)
            {
               for(int mu = -lambda; mu <= lambda; mu++)
               {
                  int dm = m - mu;

                  if(-dl <= dm && dm <= dl)
                  {
                     regularMatrix[(l * l + l + m) + (lambda * lambda + lambda + mu) * harmonicLength] =
                        regular.getHarmonic(dl * dl + dl + dm) *
                        MultipoleTranslator::localTranslationFactor(m, mu, lambda, l);
                  }
               }
            }
         }
      }
   }

   auto realToComplexMatrix = Harmonics::calcRealToComplexTransitionMatrix1D(order);
   auto complexToRealMatrix = Harmonics::calcComplexToRealTransitionMatrix1D(order);
   
   Complex alpha = make_cuComplex(1, 0);
   Complex beta = make_cuComplex(0, 0);

   std::vector<Complex> temp1(matrixElemCount);

   cblas_cgemm(
      CBLAS_ORDER::CblasRowMajor,
      CBLAS_TRANSPOSE::CblasNoTrans,
      CBLAS_TRANSPOSE::CblasNoTrans,
      harmonicLength, harmonicLength, harmonicLength,
      (float*)&alpha,
      (float*)realToComplexMatrix.data(),
      harmonicLength,
      (float*)regularMatrix.data(),
      harmonicLength,
      (float*)&beta,
      (float*)temp1.data(),
      harmonicLength);

   std::vector<Complex> temp2(matrixElemCount);
   std::vector<real> result(matrixElemCount);

   cblas_cgemm(
      CBLAS_ORDER::CblasRowMajor,
      CBLAS_TRANSPOSE::CblasNoTrans,
      CBLAS_TRANSPOSE::CblasNoTrans,
      harmonicLength, harmonicLength, harmonicLength,
      (float*)&alpha,
      (float*)temp1.data(),
      harmonicLength,
      (float*)complexToRealMatrix.data(),
      harmonicLength,
      (float*)&beta,
      (float*)temp2.data(),
      harmonicLength);

   cblas_scopy(
      matrixElemCount,
      (float*)(temp2.data()), 2,
      result.data(), 1);

   int m = harmonicLength;
   int k = harmonicLength;
   int n = 1;
   int lda = m, ldb = k, ldc = m;
   const real alpha1 = 1;
   const real beta1 = 0;
   
   std::vector<real> notTrueResult(harmonicLength);

   cblas_sgemm(CBLAS_ORDER::CblasColMajor,
               CBLAS_TRANSPOSE::CblasNoTrans,
               CBLAS_TRANSPOSE::CblasNoTrans,
               m, n, k,
               alpha1,
               result.data(), ldb, harm.data().data(), lda,
               beta1,
               notTrueResult.data(), ldc);

}

void timeForFullFMMByQuadratures()
{
   const double torusRadius = 2;
   const double torusSectionWidth = 0.2;
   Vector3 begin(3, 1, 2);
   Vector3 end(0, 0, 0);
   BasisQuadratures bq = readTetrahedronBasisQuadratures();
   int pointsCount = 1000;
   auto points = createRandomPoints({{0, 0, 0}, {2, 2, 2}}, pointsCount);

   for(size_t i = 5; i < 15; i++)
   {
      Torus torus(torusRadius, torusSectionWidth, pow(2, i + 1), 4, 4);
      auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);
      FastMultipoleSolver fmmSolver(quadratures, points, Problem::BioSavartLaplace, 1000, 100);

      auto start = std::chrono::steady_clock::now();
      fmmSolver.calcMultipoleExpansionsAtLeaves();
      fmmSolver.calclMultipoleExpansions(M2MAlg::Matrices, Device::GPU);
      auto stop = std::chrono::steady_clock::now();
      auto time1 = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      fmmSolver.calcLocalMultipoleExpansions(M2LAlg::ComplexTranslation, Device::CPU);
      stop = std::chrono::steady_clock::now();
      auto time2 = getTime(start, stop);

      std::cout << std::fixed;
      std::cout << std::setw(5) << pointsCount * quadratures.size() << " ";
      std::cout << std::scientific;
      std::cout << std::setw(8) << time1 << " " << time2 << std::endl << std::endl;
   }
}

void timeForFullFMMByPointCount()
{
   const double torusRadius = 2;
   const double torusSectionWidth = 0.2;
   Vector3 begin(3, 1, 2);
   Vector3 end(0, 0, 0);
   BasisQuadratures bq = readTetrahedronBasisQuadratures();
   Torus torus(torusRadius, torusSectionWidth, 20, 4, 4);
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

   for(size_t i = 5; i < 30; i++)
   {
      int pointsCount = pow(2, i);
      auto points = createRandomPoints({ {0, 0, 0}, {2, 2, 2} }, pointsCount);
      FastMultipoleSolver fmmSolver(quadratures, points, Problem::BioSavartLaplace, 1000, 100);

      auto start = std::chrono::steady_clock::now();
      fmmSolver.calcMultipoleExpansionsAtLeaves();
      fmmSolver.calclMultipoleExpansions(M2MAlg::Matrices, Device::GPU);
      auto stop = std::chrono::steady_clock::now();
      auto time1 = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      fmmSolver.calcLocalMultipoleExpansions(M2LAlg::ComplexTranslation, Device::CPU);
      stop = std::chrono::steady_clock::now();
      auto time2 = getTime(start, stop);

      std::cout << std::fixed;
      std::cout << std::setw(5) << pointsCount * quadratures.size() << " ";
      std::cout << std::scientific;
      std::cout << std::setw(8) << time1 << " " << time2 << std::endl;
   }
}

void timeForTallCube()
{
   const double torusRadius = 2;
   const double torusSectionWidth = 0.2;
   Vector3 begin(3, 1, 2);
   Vector3 end(0, 0, 0);
   BasisQuadratures bq = readTetrahedronBasisQuadratures();
   Torus torus(torusRadius, torusSectionWidth, 20, 4, 4);
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);
   int pointsCount = quadratures.size();

   for(size_t i = 2; i < 40; i++)
   {
      auto points = createRandomPoints({ 
         {0, 0, 0},
         {2, 2, static_cast<real>(i)}
      }, pointsCount);

      std::cout << std::scientific;

      FastMultipoleSolver fmmSolver(quadratures, points, Problem::BioSavartLaplace, 1000, 100);

      auto start = std::chrono::steady_clock::now();
      fmmSolver.calcMultipoleExpansionsAtLeaves();
      fmmSolver.calclMultipoleExpansions(M2MAlg::Matrices, Device::CPU);
      auto stop = std::chrono::steady_clock::now();
      auto time1 = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      fmmSolver.calcLocalMultipoleExpansions(M2LAlg::ComplexTranslation, Device::CPU);
      stop = std::chrono::steady_clock::now();
      auto time2 = getTime(start, stop);
         
      std::cout << std::fixed;
      std::cout << std::setw(5) << i << " ";
      std::cout << std::scientific;
      std::cout << std::setw(8) << time1 << " " << time2 << " ";

      real sumAverageError = 0;

      auto resultPoints = fmmSolver.calcB(current);

      for (auto &[point, bFmm] : resultPoints)
      {
         auto b = fmmSolver.calcB(current, point);

         sumAverageError += 100 * (bFmm - b).length() / b.length();
      }

      std::cout << std::setw(8) << sumAverageError / resultPoints.size() << std::endl;
   }
}


std::vector<ReferenceCylinderData> readCylinderData(const std::string& filename)
{
   std::vector<ReferenceCylinderData> result;

   std::ifstream fin(filename);
   std::string _;
   std::getline(fin, _);
   std::getline(fin, _);
   size_t pointId;
   
   while (fin >> pointId)
   {
      real px, py, pz;
      real bx, by, bz, bl;

      fin >> px >> py >> pz >> bx >> by >> bz >> bl;

      result.emplace_back(pointId, Vector3(px, py, pz), Vector3(bx, by, bz), bl);
   }
   
   return result;
}

void BApproximationOnCylinder()
{
   auto externalCylinderSides = readCylinderData("cylinder/ВнешнийЦилиндр.0");
   auto externalCylinderTop = readCylinderData("cylinder/ВнешнийЦилиндрВерх.0");
   auto externalCylinderBottom = readCylinderData("cylinder/ВнешнийЦилиндрНиз.0");
   auto internalCylinder = readCylinderData("cylinder/Внутренний Цилиндр.0");

   Cylinder cylinder = createCylinder();
   BasisQuadratures bq = readTriangleBasisQuadratures();

   auto BEMQuadraturesSide = math::calcBEMquadraturesFromTriangles(
      cylinder.sideTriangles(), bq, externalCylinderSides, 0);

   auto BEMQuadraturesTop = math::calcBEMquadraturesFromTriangles(
      cylinder.topTriangles(), bq, externalCylinderTop, 1);

   auto BEMQuadraturesBottom = math::calcBEMquadraturesFromTriangles(
      cylinder.bottomTriangles(), bq, externalCylinderBottom, -1);

   real sumErrorFmmSolver = 0;
   real totalTime = 0;

   for(size_t i = 0; i < internalCylinder.size(); i++)
   {
      std::cout << std::setw(5) << i << " ";

      Vector3 trueRes = internalCylinder[i].B;
      Vector3 res;

      auto start = std::chrono::steady_clock::now();
      res += math::calcBEMIntegral(internalCylinder[i].point, BEMQuadraturesSide);
      res += math::calcBEMIntegral(internalCylinder[i].point, BEMQuadraturesTop);
      res += math::calcBEMIntegral(internalCylinder[i].point, BEMQuadraturesBottom);
      auto stop = std::chrono::steady_clock::now();
      totalTime += test::getTime(start, stop);

      res.printWithWidth(std::cout, 9);
      trueRes.printWithWidth(std::cout, 9);

      real relativeError = (res - trueRes).length() / trueRes.length();

      std::cout << std::setw(10) << relativeError << std::endl;

      sumErrorFmmSolver += relativeError;
   }

   std::cout << "Average error: ";
   std::cout << std::scientific << sumErrorFmmSolver / internalCylinder.size() << std::fixed << std::endl;
   std::cout << "Total time: " << totalTime;
}

void comparisonToBEM()
{
   Cylinder cylinder = createCylinder();
   BasisQuadratures bq = readTetrahedronBasisQuadratures();

   auto externalCylinderSides = readCylinderData("cylinder/ВнешнийЦилиндр.0");
   auto externalCylinderTop = readCylinderData("cylinder/ВнешнийЦилиндрВерх.0");
   auto externalCylinderBottom = readCylinderData("cylinder/ВнешнийЦилиндрНиз.0");
   auto internalCylinder = readCylinderData("cylinder/Внутренний Цилиндр.0");

   auto BEMQuadraturesSide = math::calcBEMquadraturesFromTriangles(
      cylinder.sideTriangles(), bq, externalCylinderSides, 0);

   auto BEMQuadraturesTop = math::calcBEMquadraturesFromTriangles(
      cylinder.topTriangles(), bq, externalCylinderTop, 1);

   auto BEMQuadraturesBottom = math::calcBEMquadraturesFromTriangles(
      cylinder.bottomTriangles(), bq, externalCylinderBottom, -1);

   std::vector<BEMQuadrature> quadratures;
   quadratures.reserve(
      BEMQuadraturesSide.size() +
      BEMQuadraturesTop.size() +
      BEMQuadraturesTop.size());
      
   quadratures.insert(quadratures.end(), BEMQuadraturesSide.begin(), BEMQuadraturesSide.end());
   quadratures.insert(quadratures.end(), BEMQuadraturesTop.begin(), BEMQuadraturesTop.end());
   quadratures.insert(quadratures.end(), BEMQuadraturesBottom.begin(), BEMQuadraturesBottom.end());

   std::vector<Vector3> initialPoints(internalCylinder.size());

   for(size_t i = 0; i < internalCylinder.size(); i++)
   {
      initialPoints[i] = internalCylinder[i].point;
   }

   std::cout << "FMM SOLVER BEGIN SOLVING" << std::endl;

   auto start = std::chrono::steady_clock::now();
   FastMultipoleSolver fmmSolver(quadratures, initialPoints, Problem::BEM, 256, 128);
   fmmSolver.calclMultipoleExpansions(M2MAlg::Matrices, Device::CPU);
   fmmSolver.calcLocalMultipoleExpansions(M2LAlg::ComplexTranslation, Device::CPU);
   auto results = fmmSolver.calcBEM(current);
   auto stop = std::chrono::steady_clock::now();

   real sumErrorFmmSolver = 0;

   size_t n = internalCylinder.size();

   for(size_t i = 0; i < n; i++)
   {
      auto point = results[i].first;
      Vector3 trueRes;

      for(size_t j = 0; j < n; j++)
      {
         if(point.x == internalCylinder[j].point.x &&
            point.y == internalCylinder[j].point.y &&
            point.z == internalCylinder[j].point.z)
         {
            trueRes = internalCylinder[j].B;
            break;
         }
      }
      real errorFmm = (results[i].second - trueRes).length() / trueRes.length();

      std::cout << std::fixed << std::setw(4) << i << " ";
      
      point.printWithWidth(std::cout, 9);
      results[i].second.printWithWidth(std::cout, 9);
      
      internalCylinder[i].B.printWithWidth(std::cout, 9);

      std::cout << std::scientific << std::setw(16) << errorFmm << std::fixed;
      std::cout << std::endl;

      sumErrorFmmSolver += errorFmm;
   }

   std::cout << "Average error: ";
   std::cout << std::scientific << sumErrorFmmSolver / n << std::fixed << std::endl;
   std::cout << "Total time: " << test::getTime(start, stop);
}

int main()
{
   //NMResearch2();
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
   //layerMatrixCalculationTime(Device::CPU);
   //layerMatrixCalculationTime(Device::GPU);
   //compareWithMatrixMultiplication();

   //translationTest2();

   //FMMPrecisionTest();
   //FFMTimeTest();

   //NMResearch2();

   //timeForFullFMMByQuadratures();
   //timeForFullFMMByPointCount();
   //timeForTallCube();

   //BApproximationOnCylinder();
   comparisonToBEM();
}
