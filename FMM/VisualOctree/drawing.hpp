#pragma once
#include "GL/freeglut.h"
#include "../FMMCPU/vector3.cuh"
#include "octree.hpp"

namespace drawing
{
   struct Color
   {
      UCHAR r;
      UCHAR g;
      UCHAR b;
      UCHAR a;

      Color(UCHAR r = 0, UCHAR g = 0, UCHAR b = 0, UCHAR a = 255);
   };

   void drawPoint(Vector3 point, Color color, float size);
   void drawRectangle(Vector3 center, float halfWidth, float halfHeight, Color color);
   void drawRectangle(Vector3 a, Vector3 b, Vector3 c, Vector3 d, Color color);
   void drawOctree(const Octree& octree);
}

