#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "GL/freeglut.h"
#include "fractal.cu"

GLubyte rasters[24] = {
   0xc0, 0x00, 0xc0, 0x00, 0xc0, 0x00, 0xc0, 0x00, 0xc0, 0x00,
   0xff, 0x00, 0xff, 0x00, 0xc0, 0x00, 0xc0, 0x00, 0xc0, 0x00,
   0xff, 0xc0, 0xff, 0xc0 };

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

// Функция вывода на экран 
void Display()
{
   glClearColor(255, 255, 255, 1);
   glClear(GL_COLOR_BUFFER_BIT);

   glBitmap(10, 10, 0, 0, 10, 10, rasters);

   glFinish();
}

int main(int argc, char** argv)
{
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_RGB);
   glutInitWindowSize(DIM, DIM);
   glutCreateWindow("КГ Лабораторная работа №1");
   glutDisplayFunc(Display);
   glutReshapeFunc(Reshape);
   //glutKeyboardFunc(KeyboardLetters);
   //glutSpecialFunc(KeyboardSpecials);
   //glutMouseFunc(Mouse);

   glutMainLoop();

   return 0;
}