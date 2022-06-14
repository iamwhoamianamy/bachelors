#pragma once

//#define REAL_IS_FLOAT
#define REAL_IS_DOUBLE

#ifdef REAL_IS_FLOAT

   typedef float real;

#endif // REAL_IS_FLOAT

#ifdef REAL_IS_DOUBLE

   typedef double real;

#endif // REAL_IS_DOUBLE