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
const int IMAGE_SIZE = DIM * DIM * 3;
const int FPS = 30;

unsigned char* bitmap_pixels;

struct Info
{
   float p = CUDART_PI_F;
   float t = 0;
};

Info info;

__device__ void setPixel(unsigned char* ptr, const int offset, const unsigned char r, const unsigned char g, const unsigned char b)
{
   ptr[offset * 3 + 0] = r;
   ptr[offset * 3 + 1] = g;
   ptr[offset * 3 + 2] = b;
}

__device__ void setPixel(unsigned char* ptr, const int offset, const unsigned char value)
{
   ptr[offset * 3 + 0] = value;
   ptr[offset * 3 + 1] = value;
   ptr[offset * 3 + 2] = value;
}

__global__ void kernel(unsigned char* ptr, const float t)
{
   int x = threadIdx.x + blockIdx.x * blockDim.x;
   int y = threadIdx.y + blockIdx.y * blockDim.y;
   int offset = x + y * blockDim.x * gridDim.x;

   float fx = x - DIM / 2;
   float fy = y - DIM / 2;

   float d = sqrtf(fx * fx + fy * fy);

   unsigned char value = (unsigned char)(128.0f + 127.0f * cos(d / 10.0f - t / 7.0f) / (d / 10.0f + 1.0f));

   setPixel(ptr, offset, value);
}

void formRings()
{
   unsigned char* dev_bitmap;

   tryCudaMalloc((void**)&dev_bitmap, IMAGE_SIZE);

   dim3 blocks(DIM / 16, DIM / 16);
   dim3 threads(16, 16);

   kernel<<<blocks, threads>>>(dev_bitmap, info.t);

   tryKernelLaunch();
   tryKernelSynchronize();

   handleError(cudaMemcpy(bitmap_pixels, dev_bitmap, IMAGE_SIZE, cudaMemcpyDeviceToHost));
   cudaFree(dev_bitmap);
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

   }

   glutPostRedisplay();
}

// Функция вывода на экран 
void display()
{
   glClear(GL_COLOR_BUFFER_BIT);
   glRasterPos2i(0, 0);

   formRings();
   glDrawPixels(DIM, DIM, GL_RGB, GL_UNSIGNED_BYTE, bitmap_pixels);
   info.t += 1.0f;

   glFinish();
}

void onTimer(int millisec)
{
   glutPostRedisplay();
   glutTimerFunc(1000 / FPS, onTimer, 0);
}

void initializing()
{
   bitmap_pixels = new unsigned char[IMAGE_SIZE];
}

void exitingFunction()
{
   tryCudaLastError();
   delete[] bitmap_pixels;

   tryCudaReset();
   std::cout << "Done!";
}

int main(int argc, char** argv)
{
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_RGB);
   glutInitWindowSize(DIM, DIM);
   glutCreateWindow("Фрактал");

   glShadeModel(GL_FLAT);
   glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

   glutDisplayFunc(display);
   glutReshapeFunc(reshape);
   glutKeyboardFunc(keyboardLetters);
   atexit(exitingFunction);
   glutTimerFunc(0, onTimer, 0);

   initializing();

   glutMainLoop();

   return 0;
}