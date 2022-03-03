#include <iostream>

#include "glut_functions.hpp"
#include "window.hpp"

int main(int argc, char** argv)
{
   glut_functions::window = new Window(argc, argv, 800, 800, "Octree");
   glut_functions::window->run(30);
}
