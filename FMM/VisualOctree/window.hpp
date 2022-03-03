#pragma once
#include <iostream>
#include <vector>

#include "../FMMCPU/real.hpp"
#include "../FMMCPU/vector3.hpp"
#include "GL/freeglut.h"

class Window
{
private:
   const int FPS;
   float screenWidth;
   float screenHeight;
   std::string name;
   std::vector<Vector3> points;
public:
   Window(int argc, char** argv, int FPS,
          float screenWidth, float screenHeight, std::string name);

   void run();
   void onTimer(int millisec);
   void exitingFunction();
   void display();
   void reshape(GLint w, GLint h);

   void keyboardLetters(unsigned char key, int x, int y);
   void mouse(int button, int state, int x, int y);
   void mousePassive(int x, int y);
private:
   void initData();
};
