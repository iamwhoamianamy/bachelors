#pragma once
#include "glut_functions.hpp"

Window* glut_functions::window;

void glut_functions::registerFunctions()
{
   glutDisplayFunc(display);
   glutReshapeFunc(reshape);
   glutKeyboardFunc(keyboardLetters);
   glutMouseFunc(mouse);
   glutPassiveMotionFunc(mousePassive);
   atexit(exitingFunction);
   glutTimerFunc(0, onTimer, 0);
}

void glut_functions::onTimer(int millisec)
{
   if(window)
      window->onTimer(millisec);
}

void glut_functions::exitingFunction()
{
   if(window)
      window->exitingFunction();
}

void glut_functions::display()
{
   if(window)
      window->display();
}

void glut_functions::reshape(GLint w, GLint h)
{
   if(window)
      window->reshape(w, h);
}

void glut_functions::keyboardLetters(unsigned char key, int x, int y)
{
   if(window)
      window->keyboardLetters(key, x, y);
}

void glut_functions::mouse(int button, int state, int x, int y)
{
   if(window)
      window->mouse(button, state, x, y);
}

void glut_functions::mousePassive(int x, int y)
{
   if(window)
      window->mousePassive(x, y);
}
