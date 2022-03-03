#pragma once
#include "GL/freeglut.h"
#include "window.hpp"

namespace glut_functions
{
   extern Window* window;

   void registerFunctions();
   void onTimer(int millisec);
   void exitingFunction();
   void display();
   void reshape(GLint w, GLint h);
   void keyboardLetters(unsigned char key, int x, int y);
   void mouse(int button, int state, int x, int y);
   void mousePassive(int x, int y);
}

