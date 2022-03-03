#include "drawing.hpp"

namespace drawing
{
   Color::Color(UCHAR r, UCHAR g, UCHAR b, UCHAR a) :
      r(r), g(g), b(b), a(a) {}

   void drawPoint(Vector3 point, Color color, float size)
   {
      glColor3ub(color.r, color.g, color.b);
      glPointSize(size);
      glBegin(GL_POINTS);
      {
         glVertex3f(point.x, point.y, point.z);
      }
      glEnd();
   }

   void drawRectangle(Vector3 center, float halfWidth, float halfHeight, Color color)
   {
      glColor3ub(color.r, color.g, color.b);
      glBegin(GL_LINE_LOOP);
      {
         glVertex3f(center.x - halfWidth, center.y - halfHeight, 0);
         glVertex3f(center.x + halfWidth, center.y - halfHeight, 0);
         glVertex3f(center.x + halfWidth, center.y + halfHeight, 0);
         glVertex3f(center.x - halfWidth, center.y + halfHeight, 0);
      }
      glEnd();
   }

   void drawRectangle(Vector3 a, Vector3 b, Vector3 c, Vector3 d, Color color)
   {
      glColor3ub(color.r, color.g, color.b);
      glBegin(GL_LINE_LOOP);
      {
         glVertex3f(a.x, a.y, 0);
         glVertex3f(b.x, b.y, 0);
         glVertex3f(c.x, c.y, 0);
         glVertex3f(d.x, d.y, 0);
      }
      glEnd();
   }

   Color octreeColor(255, 255, 0);

   void drawOctree(const Octree& octree)
   {
      drawRectangle(octree.box().center,
                    octree.box().halfDimensions.x,
                    octree.box().halfDimensions.y, octreeColor);

      if(octree.isSubdivided())
      {
         for(auto child : octree.children())
         {
            drawOctree(*child);
         }
      }
   }
}
