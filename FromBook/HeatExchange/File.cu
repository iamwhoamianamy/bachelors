#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "GL/freeglut.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cudahelper.cuh"
#include "math_constants.h"

const int DIM = 1024;
const int GRID_SIZE = DIM * DIM;
const int IMAGE_SIZE = DIM * DIM * 3;
const int FPS = 30;
__device__ const float SPEED = 0.25f;
const float MAX_TEMP = 1.0f;
const float MIN_TEMP = 0.0001f;

const dim3 blocks(DIM / 16, DIM / 16);
const dim3 threads(16, 16);

struct DataBlock
{
   unsigned char* output_pixels;
   unsigned char* dev_pixels;
   float* dev_inSrc;
   float* dev_outSrc;
   float* dev_constSrc;
   float* constSrc;
};

DataBlock data;

__device__ void setPixel(unsigned char* ptr, int offset, const unsigned char r, const unsigned char g, const unsigned char b)
{
   offset *= 3;
   ptr[offset + 0] = r;
   ptr[offset + 1] = g;
   ptr[offset + 2] = b;
}

__device__ void setPixel(unsigned char* ptr, int offset, const unsigned char value)
{
   offset *= 3;
   ptr[offset + 0] = value;
   ptr[offset + 1] = value;
   ptr[offset + 2] = value;
}

__global__ void clearArray(float* arr)
{
   int x = threadIdx.x + blockIdx.x * blockDim.x;
   int y = threadIdx.y + blockIdx.y * blockDim.y;
   int offset = x + y * blockDim.x * gridDim.x;

   arr[offset] = 0;
}

__global__ void copyHeatersKernel(float* inPtr, const float* heatersMap)
{
   int x = threadIdx.x + blockIdx.x * blockDim.x;
   int y = threadIdx.y + blockIdx.y * blockDim.y;
   int offset = x + y * blockDim.x * gridDim.x;

   if(heatersMap[offset] != 0)
      inPtr[offset] = heatersMap[offset];
}

__global__ void blendKernel(float* outSrc, const float* inSrc)
{
   int x = threadIdx.x + blockIdx.x * blockDim.x;
   int y = threadIdx.y + blockIdx.y * blockDim.y;
   int offset = x + y * blockDim.x * gridDim.x;

   int left = x > 0 ? offset - 1 : offset;
   int right = x + 1 < DIM ? offset + 1 : offset;

   int bot = y > 0 ? offset - DIM : offset;
   int top = y + 1 < DIM ? offset + DIM : offset;

   outSrc[offset] = inSrc[offset] + SPEED * (inSrc[top] + inSrc[bot] + inSrc[left] + inSrc[right] - 4 * inSrc[offset]);
}

__global__ void floatToColor(const float* values, unsigned char* colors)
{
   int x = threadIdx.x + blockIdx.x * blockDim.x;
   int y = threadIdx.y + blockIdx.y * blockDim.y;
   int offset = x + y * blockDim.x * gridDim.x;

   setPixel(colors, offset, (unsigned char)(values[offset] * 255));
}

handler freeMemory();

void initMemory()
{
   tryCudaMalloc((void**)&data.dev_pixels, IMAGE_SIZE * sizeof(float));
   tryCudaMalloc((void**)&data.dev_inSrc, GRID_SIZE * sizeof(float));
   tryCudaMalloc((void**)&data.dev_outSrc, GRID_SIZE * sizeof(float));
   tryCudaMalloc((void**)&data.dev_constSrc, GRID_SIZE * sizeof(float));

   data.output_pixels = new unsigned char[IMAGE_SIZE];
   data.constSrc = new float[GRID_SIZE];
}

void formHeatersMap()
{
   for(size_t i = 0; i < GRID_SIZE; i++)
   {
      data.constSrc[i] = 0;
   }

   /*data.constSrc[100 + DIM * 100] = MAX_TEMP;
   data.constSrc[110 + DIM * 100] = MAX_TEMP;
   data.constSrc[120 + DIM * 100] = MAX_TEMP;
   data.constSrc[130 + DIM * 100] = MAX_TEMP;
   data.constSrc[140 + DIM * 100] = MAX_TEMP;
   data.constSrc[300 + DIM * 300] = MAX_TEMP;*/

   handleError(cudaMemcpy(data.dev_constSrc, data.constSrc, GRID_SIZE, cudaMemcpyHostToDevice));
}

const int REPEATS = 10;

void formHeatmap()
{
   

   for(size_t i = 0; i < REPEATS; i++)
   {
      copyHeatersKernel << <blocks, threads >> > (data.dev_inSrc, data.dev_constSrc);
      tryKernelLaunch();
      tryKernelSynchronize();

      blendKernel << <blocks, threads >> > (data.dev_outSrc, data.dev_inSrc);
      tryKernelLaunch();
      tryKernelSynchronize();

      std::swap(data.dev_inSrc, data.dev_outSrc);
   }

   floatToColor << <blocks, threads >> > (data.dev_inSrc, data.dev_pixels);
   tryKernelLaunch();
   tryKernelSynchronize();

   handleError(cudaMemcpy(data.output_pixels, data.dev_pixels, IMAGE_SIZE * sizeof(unsigned char), cudaMemcpyDeviceToHost));
}

handler freeMemory()
{
   cudaFree(data.dev_pixels);
   cudaFree(data.dev_inSrc);
   cudaFree(data.dev_outSrc);
   cudaFree(data.dev_constSrc);

   delete[] data.output_pixels;
   delete[] data.constSrc;
}

// Функция изменения размеров окна
void reshape(GLint w, GLint h)
{
   glViewport(0, 0, DIM, DIM);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glOrtho(0, DIM, 0, DIM, -1.0, 1.0);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
}

// Функция обработки сообщений от клавиатуры 1
void keyboardLetters(unsigned char key, int x, int y)
{
   switch(key)
   {
      case 'r':
      {
         for(size_t i = 0; i < GRID_SIZE; i++)
         {
            data.constSrc[i] = 0;
         }

         handleError(cudaMemcpy(data.dev_constSrc, data.constSrc, GRID_SIZE * sizeof(float), cudaMemcpyHostToDevice));
         break;
      }
      case 'c':
      {
         clearArray << <blocks, threads >> > (data.dev_inSrc);
         tryKernelLaunch();
         tryKernelSynchronize();

         break;
      }
   }

   //glutPostRedisplay();
}

// Функция обработки сообщения от мыши
void Mouse(int button, int state, int x, int y)
{
   // Клавиша была нажата, но не отпущена
   if(state != GLUT_DOWN) return;

   // Новая точка по левому клику
   if(button == GLUT_LEFT_BUTTON)
   {
      data.constSrc[x + DIM * (DIM - y)] = MAX_TEMP;

      //cudaFree(data.dev_constSrc);
      //tryCudaMalloc((void**)&data.dev_constSrc, GRID_SIZE * sizeof(float));
      handleError(cudaMemcpy(data.dev_constSrc, data.constSrc, GRID_SIZE * sizeof(float), cudaMemcpyHostToDevice));
   }

   // Новая точка по левому клику
   if(button == GLUT_RIGHT_BUTTON)
   {
      data.constSrc[x + DIM * (DIM - y)] = MIN_TEMP;

      //cudaFree(data.dev_constSrc);
      //tryCudaMalloc((void**)&data.dev_constSrc, GRID_SIZE * sizeof(float));
      handleError(cudaMemcpy(data.dev_constSrc, data.constSrc, GRID_SIZE * sizeof(float), cudaMemcpyHostToDevice));
   }
}

// Функция вывода на экран 
void display()
{
   glClear(GL_COLOR_BUFFER_BIT);
   glRasterPos2i(0, 0);

   formHeatmap();
   glDrawPixels(DIM, DIM, GL_RGB, GL_UNSIGNED_BYTE, data.output_pixels);

   glFinish();
}

void onTimer(int millisec)
{
   glutPostRedisplay();
   glutTimerFunc(1000 / FPS, onTimer, 0);
}

void exitingFunction()
{
   freeMemory();
   //tryCudaLastError();
   tryCudaReset();
   std::cout << "Done!";
}

int main(int argc, char** argv)
{
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_RGB);
   glutInitWindowSize(DIM, DIM);
   glutCreateWindow("Нагрев");

   glShadeModel(GL_FLAT);
   glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

   glutDisplayFunc(display);
   glutReshapeFunc(reshape);
   glutKeyboardFunc(keyboardLetters);
   glutMouseFunc(Mouse);
   atexit(exitingFunction);
   glutTimerFunc(0, onTimer, 0);

   initMemory();
   glutMainLoop();

   return 0;
}