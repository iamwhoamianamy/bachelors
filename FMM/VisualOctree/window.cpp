#include "window.hpp"
#include "drawing.hpp"
#include "glut_functions.hpp"
#include "../FMMCPU/vector3.cpp"

namespace glf = glut_functions;

Window::Window(int argc, char** argv, float screenWidth, float screenHeight, std::string name) :
   screenWidth(screenWidth), screenHeight(screenHeight), name(name)
{
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_RGB);
   glutInitWindowSize(screenWidth, screenHeight);
   glutCreateWindow(name.c_str());
      
   glf::registerFunctions();
   initData();
}

void Window::initData()
{
   points = std::vector<Vector3>();
   auto vec = Vector3(200, 200);
   points.push_back(vec);
}

void Window::run(int FPS)
{
   this->FPS = FPS;
   glutMainLoop();
}

void Window::onTimer(int millisec)
{
   glutPostRedisplay();
   glutTimerFunc(1000 / FPS, glut_functions::onTimer, 0);
}

void Window::exitingFunction()
{

   std::cout << "Done!";
}

void Window::reshape(GLint w, GLint h)
{
   screenWidth = w;
   screenHeight = h;

   glViewport(0, 0, screenWidth, screenHeight);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glOrtho(0, screenWidth, screenHeight, 0, -1.0, 1.0);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
}

void Window::keyboardLetters(unsigned char key, int x, int y)
{

}

void Window::mouse(int button, int state, int x, int y)
{
   points.push_back(Vector3(x, y, 0));
}

void Window::mousePassive(int x, int y)
{
   mousePos.x = x;
   mousePos.y = y;
}

void Window::display()
{
   glClearColor(0, 0, 0, 255);
   glClear(GL_COLOR_BUFFER_BIT);

   auto pointColor = drawing::Color(255, 255, 255);

   for(auto& point : points)
   {
      drawing::DrawPoint(point, pointColor, 20);
   }

   drawing::DrawRectangle(mousePos, 20, 20, pointColor);

   glFinish();
}