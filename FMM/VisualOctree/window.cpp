#include "window.hpp"
#include "glut_functions.hpp"

namespace glf = glut_functions;

Window::Window(int argc, char** argv, int FPS, float screenWidth, float screenHeight, std::string name) :
   FPS(FPS), screenWidth(screenWidth), screenHeight(screenHeight), name(name)
{
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_RGB);
   glutInitWindowSize(screenWidth, screenHeight);
   glutCreateWindow(name.c_str());

   glShadeModel(GL_FLAT);
   glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
   
   glf::registerFunctions();
   glutMainLoop();
}

void Window::OnTimer(int millisec)
{
   glutPostRedisplay();
   glutTimerFunc(1000 / FPS, glut_functions::onTimer, 0);
}

void Window::ExitingFunction()
{

   std::cout << "Done!";
}

void Window::Reshape(GLint w, GLint h)
{
   glViewport(0, 0, screenWidth, screenHeight);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glOrtho(0, screenWidth, 0, screenHeight, -1.0, 1.0);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
}

void Window::KeyboardLetters(unsigned char key, int x, int y)
{

}

void Window::Mouse(int button, int state, int x, int y)
{

}

void Window::MousePassive(int x, int y)
{

}


void Window::Display()
{
   glClearColor(100, 0, 0, 255);
   glClear(GL_COLOR_BUFFER_BIT);



   glFinish();
}