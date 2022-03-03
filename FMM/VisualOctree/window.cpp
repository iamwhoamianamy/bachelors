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
   if(state == 0)
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

   glEnable(GL_POINT_SMOOTH);
      
   for(auto& point : points)
   {
      drawing::drawPoint(point, pointColor, 5);
   }

   Box mouseBox(mousePos, Vector3(40, 40, 40));
   drawing::drawRectangle(mousePos, mouseBox.halfDimensions.x, mouseBox.halfDimensions.y, pointColor);

   Octree octree(Box(Vector3(screenWidth / 2, screenHeight / 2, 0),
                 Vector3(screenWidth / 2, screenHeight / 2, 50)), 1);
   std::vector<Vector3*> foundPoints;

   octree.insert(points);
   octree.quarry(mouseBox, foundPoints);

   auto foundColor = drawing::Color(255, 0, 0);

   for(auto point : foundPoints)
   {
      drawing::drawPoint(*point, foundColor, 5);
   }

   drawing::drawOctree(octree);

   glFinish();
}