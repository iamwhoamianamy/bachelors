#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "GL/freeglut.h"
#include "bitmap_image.hpp"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cudahelper.h"
#include "math_constants.h"

const int DIM = 1000;
const int IMAGE_SIZE = DIM * DIM * 3;
__device__ const float SCALE = 1.5;

unsigned char *bitmap_pixels;

struct FractalInfo
{
   float p = CUDART_PI_F;
};

FractalInfo fractalInfo;

__device__ const float r = 0.7885;
__device__ const float i = 0.156;


struct Complex
{
   float r;
   float i;

   __device__ Complex(float a, float b) : r(a), i(b) {}

   __device__ float MagSquared()
   {
      return r * r + i * i;
   }

   __device__ Complex operator * (const Complex& rhs)
   {
      return Complex(r * rhs.r - i * rhs.i, i * rhs.r + r * rhs.i);
   }

   __device__ Complex operator + (const Complex& rhs)
   {
      return Complex(r + rhs.r, i + rhs.i);
   }
};

__device__ void SetPixel(unsigned char* ptr, const int offset, const unsigned char r, const unsigned char g, const unsigned char b)
{
   ptr[offset * 3 + 0] = r;
   ptr[offset * 3 + 1] = g;
   ptr[offset * 3 + 2] = b;
}

/// <returns>Value from 0 to 255</returns>
__device__ int Julia(const int x, const int y, const float p)
{
   float jx = SCALE * (float)(DIM / 2 - x) / (DIM / 2);
   float jy = SCALE * (float)(DIM / 2 - y) / (DIM / 2);

   Complex c(r * cos(p), r * sin(p));
   Complex a(jx, jy);

   size_t i;
   const int max_iter = 255;
   for(i = 1; i < max_iter; i++)
   {
      a = a * a + c;
      if(a.MagSquared() > 1000)
         break;
   }

   return max_iter / i;
}

__device__ struct RGBStruct
{
   float R;
   float G;
   float B;

   __device__ RGBStruct() : R(0), G(0), B(0) {}
   __device__ RGBStruct(float r, float g, float b) : R(r), G(g), B(b) {}
};

__device__ RGBStruct HSVToRGB(float h, float s, float v)
{
   float c = 0.0f, m = 0.0f, x = 0.0f;
   RGBStruct color;
   
   c = v * s;
   x = c * (1.0 - abs((int)(h / 60.0) % 2 - 1.0));
   m = v - c;

   if(h >= 0.0 && h < 60.0)
   {
      color = RGBStruct(c + m, x + m, m);
   }
   else if(h >= 60.0 && h < 120.0)
   {
      color = RGBStruct(x + m, c + m, m);
   }
   else if(h >= 120.0 && h < 180.0)
   {
      color = RGBStruct(m, c + m, x + m);
   }
   else if(h >= 180.0 && h < 240.0)
   {
      color = RGBStruct(m, x + m, c + m);
   }
   else if(h >= 240.0 && h < 300.0)
   {
      color = RGBStruct(x + m, m, c + m);
   }
   else if(h >= 300.0 && h < 360.0)
   {
      color = RGBStruct(c + m, m, x + m);
   }
   else
   {
      color = RGBStruct(m, m, m);
   }

   return color;
}

__global__ void kernel(unsigned char* ptr, const float p)
{
   int x = blockIdx.x;
   int y = blockIdx.y;
   int offset = x + y * gridDim.x;

   float julia_value = Julia(x, y, p);

   //RGBStruct color = HSVToRGB(julia_value * 4.0f, 1, 1);
   //SetPixel(ptr, offset, (unsigned char)(color.R) * 255, (unsigned char)(color.G) * 255, (unsigned char)(color.B) * 255);
   SetPixel(ptr, offset, julia_value, 0, 0);
}

void FormFractal()
{
   unsigned char* dev_bitmap;

   TryCudaMalloc((void**)&dev_bitmap, IMAGE_SIZE);

   dim3 grid(DIM, DIM);

   kernel << <grid, 1 >> > (dev_bitmap, fractalInfo.p);

   HANDLE_ERROR(cudaMemcpy(bitmap_pixels, dev_bitmap, IMAGE_SIZE, cudaMemcpyDeviceToHost));
   cudaFree(dev_bitmap);
}

// Функция изменения размеров окна
void Reshape(GLint w, GLint h)
{
   glViewport(0, 0, DIM, DIM);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glOrtho(0, DIM, 0, DIM, -1.0, 1.0);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
}

// Функция обработки сообщений от клавиатуры 1
void KeyboardLetters(unsigned char key, int x, int y)
{
   switch(key)
   {
      //case '+': fractalInfo.scale += 0.1; FormFractal(); MakeCheckImage(); break;
      //case '-': fractalInfo.scale -= 0.1; FormFractal(); MakeCheckImage(); break;
      case ',':
      {
         fractalInfo.p += 0.01f;

         if(fractalInfo.p >= CUDART_PI_F * 2.0f)
            fractalInfo.p = 0;

         FormFractal();
         break;
      }
      case '.':
      {
         fractalInfo.p -= 0.01f;

         if(fractalInfo.p <= 0)
            fractalInfo.p = CUDART_PI_F * 2.0f;

         FormFractal();
         break;
      }
   }

   glutPostRedisplay();
}

// Функция вывода на экран 
void Display()
{
   glClear(GL_COLOR_BUFFER_BIT);
   glRasterPos2i(0, 0);

   glDrawPixels(DIM, DIM, GL_RGB, GL_UNSIGNED_BYTE, bitmap_pixels);

   glFinish();
}

int main(int argc, char** argv)
{
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_RGB);
   glutInitWindowSize(DIM, DIM);
   glutCreateWindow("Фрактал");

   glShadeModel(GL_FLAT);
   glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

   bitmap_pixels = new unsigned char[IMAGE_SIZE];
   FormFractal();

   glutDisplayFunc(Display);
   glutReshapeFunc(Reshape);
   glutKeyboardFunc(KeyboardLetters);
   //glutSpecialFunc(KeyboardSpecials);
   //glutMouseFunc(Mouse);

   glutMainLoop();

   delete[] bitmap_pixels;
   return 0;
}