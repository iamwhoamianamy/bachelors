#pragma once

enum class M2MAlg
{
   NoTranslation = 0,
   ComplexTranslation,
   RealTranslation,
   Layers,
   Matrices,
};

enum class Device
{
   CPU = 0,
   GPU,
   Adaptive
};


enum class M2LAlg
{
   ComplexTranslation = 0,
   Matrices,
};

enum class Problem
{
   BioSavartLaplace = 0,
   BEM,
};