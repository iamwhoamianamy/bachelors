#include "window.hpp"
#include "glut_functions.hpp"
#include "../FMMCPU/vector3.cpp"

namespace glf = glut_functions;

Window::Window(int argc, char** argv, int FPS, float screenWidth, float screenHeight, std::string name) :
   FPS(FPS), screenWidth(screenWidth), screenHeight(screenHeight), name(name)
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

void Window::run()
{
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

}


void Window::display()
{
   glClearColor(0, 0, 0, 255);
   glClear(GL_COLOR_BUFFER_BIT);

   glColor3ub(255, 255, 255);
   glPointSize(10);
   for (auto &point : points)
   {
      glBegin(GL_POINTS);
      {
         glVertex2f(point.x, point.y);
      }
      glEnd();
   }

   glFinish();
}