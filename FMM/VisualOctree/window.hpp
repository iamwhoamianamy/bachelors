#pragma once
#include <iostream>
#include <vector>

#include "GL/freeglut.h"

class Window
{
private:
   const int FPS;
   float screenWidth;
   float screenHeight;
   std::string name;
public:
   Window(int argc, char** argv, int FPS,
          float screenWidth, float screenHeight, std::string name);

   void OnTimer(int millisec);
   void ExitingFunction();
   void Display();
   void Reshape(GLint w, GLint h);

   void KeyboardLetters(unsigned char key, int x, int y);
   void Mouse(int button, int state, int x, int y);
   void MousePassive(int x, int y);
};
