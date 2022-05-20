#pragma once

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
   ComplexTranslation = 0,
   Matrices,
};