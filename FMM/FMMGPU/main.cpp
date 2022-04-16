#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>

#include "vector3.cuh"
#include "math.hpp"
#include "integration.hpp"
#include "harmonics.hpp"
#include "multipole_solver.hpp"
#include "translation_algorithms.hpp"
#include "matrix_mult.hpp"
#include "testing_helpers.hpp"

using namespace math;
using namespace test;
const real current = 5;

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
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

   auto telmaResults = readTelmaResults("results/telma_results.txt");
   real sumError = 0;
   size_t n = telmaResults.size();

   for(size_t i = 0; i < n; i++)
   {
      auto point = telmaResults[i].first;
      auto telmaB = telmaResults[i].second * math::mu0;

      Vector3 myB = math::calcBioSavartLaplace(current, point, quadratures);
      real error = 100 * (myB - telmaB).length() / telmaB.length();

      std::cout << std::fixed << std::setw(3) << i << " ";
      point.printWithWidth(std::cout, 6);
      std::cout << std::setw(16) << error << std::endl;

      sumError += error;
   }

   std::cout << sumError / n << std::endl;
}

void comparisonToTelmaWithoutTranslation()
{
   Torus torus = createTorus();
   BasisQuadratures bq = readBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

   auto telmaResults = readTelmaResults("results/telma_results.txt");
   MultipoleSolver multipoleSolver(quadratures);
   multipoleSolver.calcLocalMultipolesWithoutTranslation();
   real sumError = 0;
   size_t n = telmaResults.size();

   for(size_t i = 0; i < n; i++)
   {
      auto point = telmaResults[i].first;
      auto telmaB = telmaResults[i].second * math::mu0;

      Vector3 myB = multipoleSolver.calcB(current, point);
      real error = 100 * (myB - telmaB).length() / telmaB.length();

      std::cout << std::fixed << std::setw(3) << i << " ";
      point.printWithWidth(std::cout, 6);
      std::cout << std::setw(16) << error << std::endl;

      sumError += error;
   }

   std::cout << sumError / n << std::endl;
}

void comparisonToTelmaWithTranslation()
{
   Torus torus = createTorus();
   BasisQuadratures bq = readBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

   auto telmaResults = readTelmaResults("results/telma_results.txt");
   MultipoleSolver multipoleSolver(quadratures);
   multipoleSolver.calcLocalMultipolesWithRealTranslation();
   real sumError = 0;
   size_t n = telmaResults.size();

   for(size_t i = 0; i < n; i++)
   {
      auto point = telmaResults[i].first;
      auto telmaB = telmaResults[i].second * math::mu0;

      Vector3 myB = multipoleSolver.calcB(current, point);
      real error = 100 * (myB - telmaB).length() / telmaB.length();

      std::cout << std::fixed << std::setw(3) << i << " ";
      point.printWithWidth(std::cout, 6);
      std::cout << std::setw(16) << error << std::endl;

      sumError += error;
   }

   std::cout << sumError / n << std::endl;
}

void comparisonBetweenMethodsOnPrecision()
{
   Torus torus = createTorus();
   BasisQuadratures bq = readBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);
   MultipoleSolver multipoleSolver(quadratures);

   Vector3 point(3, 1, 2);

   Vector3 byIntegration = math::calcBioSavartLaplace(current, point, quadratures);
   
   multipoleSolver.calcLocalMultipolesWithoutTranslation();
   Vector3 byMultipolesWithoutTranslation = multipoleSolver.calcB(current, point);

   multipoleSolver.calcLocalMultipolesWithComplexTranslation();
   Vector3 byMultipolesWithComplexTranslation = multipoleSolver.calcB(current, point);

   multipoleSolver.calcLocalMultipolesWithRealTranslation();
   Vector3 byMultipolesWithRealTranslation = multipoleSolver.calcB(current, point);

   multipoleSolver.calcLocalMultipolesWithLayers(false);
   Vector3 byMultipolesWithLayersCPU = multipoleSolver.calcB(current, point);

   multipoleSolver.calcLocalMultipolesWithLayers(true);
   Vector3 byMultipolesWithLayersGPU = multipoleSolver.calcB(current, point);

   multipoleSolver.calcLocalMultipolesWithMatrices(false);
   Vector3 byMultipolesWithMatricesCPU = multipoleSolver.calcB(current, point);

   multipoleSolver.calcLocalMultipolesWithMatrices(true);
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

void translationTest()
{
   Vector3 point1(3, 1, 2);
   auto r1 = Harmonics::calcRegularSolidHarmonics(10, point1);
   auto c1 = Harmonics::realToComplex(r1);

   Vector3 point2(4, 5, 1);
   auto r2 = Harmonics::calcRegularSolidHarmonics(10, point2);
   auto c2 = Harmonics::realToComplex(r2);

   auto t =  Harmonics::complexToReal(Harmonics::translate(c1, c2));
}

void calculationTimeForLocalMultipoles()
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
   std::cout << std::endl;

   std::cout << std::fixed;

   for(size_t i = 1; i < 3; i++)
   {
      int octreeLeafCapacity = pow(2, i);
      MultipoleSolver multipoleSolverWithoutT(quadratures, octreeLeafCapacity);
      MultipoleSolver multipoleSolverWithComplexT(quadratures, octreeLeafCapacity);
      MultipoleSolver multipoleSolverWithRealT(quadratures, octreeLeafCapacity);
      MultipoleSolver multipoleSolverWithLayersCPU(quadratures, octreeLeafCapacity);
      MultipoleSolver multipoleSolverWithLayersGPU(quadratures, octreeLeafCapacity);

      auto start = std::chrono::steady_clock::now();
      //multipoleSolverWithoutT.calcLocalMultipolesWithoutTranslation();
      auto stop = std::chrono::steady_clock::now();
      double timeWithoutTranslation = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      //multipoleSolverWithComplexT.calcLocalMultipolesWithComplexTranslation();
      stop = std::chrono::steady_clock::now();
      double timeWithComplexTranslation = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      //multipoleSolverWithRealT.calcLocalMultipolesWithRealTranslation();
      stop = std::chrono::steady_clock::now();
      double timeWithRealTranslation = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolverWithLayersCPU.calcLocalMultipolesWithLayers(false);
      stop = std::chrono::steady_clock::now();
      double timeWithLayersCPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolverWithLayersGPU.calcLocalMultipolesWithLayers(true);
      stop = std::chrono::steady_clock::now();
      double timeWithLayersGPU = getTime(start, stop);

      std::cout << std::setw(w) << octreeLeafCapacity;
      std::cout << std::setw(w) << timeWithoutTranslation;
      std::cout << std::setw(w) << timeWithComplexTranslation;
      std::cout << std::setw(w) << timeWithRealTranslation;
      std::cout << std::setw(w) << timeWithLayersCPU;
      std::cout << std::setw(w) << timeWithLayersGPU << std::endl;
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

std::pair<double, Vector3>  timeForMultipolesWithoutTranslation(
   const std::vector<Vector3>& points,
   std::vector<Quadrature>& quadratures)
{
   Vector3 res;
   MultipoleSolver multipoleSolver(quadratures);
   multipoleSolver.calcLocalMultipolesWithoutTranslation();

   auto start = std::chrono::steady_clock::now();

   for(size_t p = 0; p < points.size(); p++)
   {
      res += multipoleSolver.calcB(current, points[p]);
   }

   auto stop = std::chrono::steady_clock::now();

   return std::pair<double, Vector3>
      (std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count() * 1e-9, res);
}

std::pair<double, Vector3>  timeForMultipolesWithTranslation(
   const std::vector<Vector3>& points,
   std::vector<Quadrature>& quadratures)
{
   Vector3 res;
   MultipoleSolver multipoleSolver(quadratures);
   multipoleSolver.calcLocalMultipolesWithRealTranslation();

   auto start = std::chrono::steady_clock::now();

   for(size_t p = 0; p < points.size(); p++)
   {
      res += multipoleSolver.calcB(current, points[p]);
   }

   auto stop = std::chrono::steady_clock::now();

   return std::pair<double, Vector3>
      (std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count() * 1e-9, res);
}

void timeResearchForMorePoints()
{
   Torus torus = createTorus();
   BasisQuadratures bq = readBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);
   
   Vector3 begin(3, 1, 2);
   Vector3 end(0, 0, 0);
   
   for(size_t i = 2; i < 14; i++)
   {
      int pointsCount = pow(2, i);
      auto points = createPoints(begin, end, pointsCount);
      
      std::cout << std::fixed;
      std::cout << std::setw(5) << pointsCount << " ";
      std::cout << std::scientific;
      std::cout << std::setw(8) << wholeTimeForIntegrals(points, quadratures).first << " ";
      std::cout << std::setw(8) << timeForMultipolesWithoutTranslation(points, quadratures).first << " ";
      std::cout << std::setw(8) << timeForMultipolesWithTranslation(points, quadratures).first << std::endl;
   }
}

void NMResearch()
{
   const double torusRadius = 2;
   const double torusSectionWidth = 0.2;
   Vector3 begin(3, 1, 2);
   Vector3 end(0, 0, 0);
   BasisQuadratures bq = readBasisQuadratures();

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
      std::cout << std::setw(8) << timeForMultipolesWithoutTranslation(points, quadratures).first << " ";
      std::cout << std::setw(8) << timeForMultipolesWithTranslation(points, quadratures).first << std::endl;
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

   multipoleSolverCPU.calcLocalMultipolesWithLayers(false);
   Vector3 byMultipolesWithLayersCPU = multipoleSolverCPU.calcB(current, point);

   multipoleSolverGPU.calcLocalMultipolesWithLayers(true);
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
      multipoleSolverCPU.calcLocalMultipolesWithLayers(0);
      auto stop = std::chrono::steady_clock::now();
      double timeWithCPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolverGPU.calcLocalMultipolesWithLayers(1);
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
      //MultipoleSolver multipoleSolverCPU(quadratures, octreeLeafCapacity);
      MultipoleSolver multipoleSolverGPU(quadratures, octreeLeafCapacity);

      std::cout << std::setw(w) << "leaf capacity:";
      std::cout << std::setw(w) << octreeLeafCapacity << std::endl;

      auto start = std::chrono::steady_clock::now();
      //multipoleSolverCPU.calcLocalMultipolesWithMatrices(0);
      auto stop = std::chrono::steady_clock::now();
      double timeWithCPU = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolverGPU.calcLocalMultipolesWithMatrices(1);
      stop = std::chrono::steady_clock::now();
      double timeWithGPU = getTime(start, stop);

      //std::cout << std::setw(w) << "total time CPU:";
      //std::cout << std::setw(w) << timeWithCPU << std::endl;
      std::cout << std::setw(w) << "total time GPU:";
      std::cout << std::setw(w) << timeWithGPU << std::endl;
   }
}

void layerMatrixCalculationTime(bool useGPU = 0)
{
   Torus torus = createTorus();
   BasisQuadratures bq = readBasisQuadratures();
   auto quadratures = math::tetrahedraToQuadratures(torus.tetrahedra, bq);

   Vector3 point(3, 1, 2);

   size_t w = 15;

   for(size_t i = 1; i < 2; i++)
   {
      int octreeLeafCapacity = pow(2, i);
      MultipoleSolver multipoleSolverLayers(quadratures, octreeLeafCapacity);
      MultipoleSolver multipoleSolverMatrices(quadratures, octreeLeafCapacity);

      std::cout << std::setw(w) << "leaf capacity:";
      std::cout << std::setw(w) << octreeLeafCapacity << std::endl;

      auto start = std::chrono::steady_clock::now();
      multipoleSolverLayers.calcLocalMultipolesWithLayers(useGPU);
      auto stop = std::chrono::steady_clock::now();
      double timeWithLayers = getTime(start, stop);

      start = std::chrono::steady_clock::now();
      multipoleSolverMatrices.calcLocalMultipolesWithMatrices(useGPU);
      stop = std::chrono::steady_clock::now();
      double timeWithMatrices = getTime(start, stop);

      
      std::cout << "device" << std::setw(w) << (useGPU ? "GPU" : "CPU") << std::endl;
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
   //MultipoleSolver multipoleSolverCPU(quadratures);
   MultipoleSolver multipoleSolverGPU(quadratures);

   Vector3 point(3, 1, 2);

   Vector3 byIntegration = math::calcBioSavartLaplace(current, point, quadratures);

   //multipoleSolverCPU.calcLocalMultipolesWithMatrices(false);
   //Vector3 byMultipolesWithLayersCPU = multipoleSolverCPU.calcB(current, point);

   multipoleSolverGPU.calcLocalMultipolesWithMatrices(true);
   Vector3 byMultipolesWithLayersGPU = multipoleSolverGPU.calcB(current, point);

   std::cout << std::setw(20) << "point ";
   point.printWithWidth(std::cout, 6);
   std::cout << std::scientific << std::endl;
   std::cout << std::setw(40) << "integration " << byIntegration << std::endl;
   //std::cout << std::setw(40) << "multipoles with matrices CPU" << byMultipolesWithLayersCPU << std::endl;
   std::cout << std::setw(40) << "multipoles with matrices GPU" << byMultipolesWithLayersGPU << std::endl;
}

int main()
{
   //NMResearch();
   //timeResearchForMorePoints();
   //comparisonToTelmaWithTranslation();
   comparisonBetweenMethodsOnPrecision();
   //translationTest();
   //calculationTimeForLocalMultipoles();
   //layerCalculationsPrecision();
   //matrixCalculationsPrecision();

   //layerCalculationTime();
   //matrixCalculationTime();
   //layerMatrixCalculationTime(0);
   //layerMatrixCalculationTime(1);

   //compareWithMatrixMultiplication();

   /*Vector3 point(3, 1, 2);
   auto a = Harmonics::calcRegularSolidHarmonics(10, point);
   auto b = Harmonics::calcRegularSolidHarmonics(10, point);

   auto ac = Harmonics::realToComplex(a);
   auto bc = Harmonics::realToComplex(b);

   Harmonics::mult(ac, 1.0 / Complex(0, 1));
   Harmonics::multToIDiv(bc);*/

}
