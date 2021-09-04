#pragma once
#include <vector>

using namespace std;

class NotImplementedExeption{};

class QuadPoints
{
public:
   vector<double> x;
   vector<double> y;
   vector<double> w;

   int order;

   QuadPoints(int order) : order(order)
   {
      switch(order)
      {
         case 6:
         {
            x = { 0.659027622374092, 0.659027622374092, 0.231933368553031, 0.231933368553031, 0.109039009072877, 0.109039009072877 };
            y = { 0.231933368553031, 0.109039009072877, 0.659027622374092, 0.109039009072877, 0.659027622374092, 0.231933368553031 };
            w = { 0.166666666666667, 0.166666666666667, 0.166666666666667, 0.166666666666667, 0.166666666666667, 0.166666666666667 };
            break;
         }
         default:
         {
            throw NotImplementedExeption();
         }
      }
   }
};

struct Vertix
{
   double x;
   double y;

   Vertix(double x, double y) : x(x), y(y) {}
};

class Triangle
{
private:
   Vertix a;
   Vertix b;
   Vertix c;

   double area;

   Triangle();

   void SetArea()
   {
      area = 0.5 * abs((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
   }

public:
   Triangle(double ax, double ay,
            double bx, double by,
            double cx, double cy) : a(ax, ay), b(bx, by), c(cx, cy)
   {
      SetArea();
   }

   Triangle(Vertix a, Vertix b, Vertix c) : a(a), b(b), c(c)
   {
      SetArea();
   }

   double GetArea() const { return area; }

   vector<double> Xs() const
   {
      return { a.x, b.x, c.x };
   }

   vector<double> Ys() const
   {
      return { a.y, b.y, c.y };
   }

   class RangeExeption{};

   Vertix operator[](int i) const
   {
      switch(i)
      {
         case 0: return a;
         case 1: return b;
         case 2: return c;
         default: throw RangeExeption();
      }
   }

   double Xst(double s, double t) const
   {
      return a.x + (b.x - a.x) * s + (c.x - a.x) * t;
   }

   double Yst(double s, double t) const
   {
      return a.y + (b.y - a.y) * s + (c.y - a.y) * t;
   }

   Vertix Vst(double s, double t) const
   {
      return Vertix(Xst(s, t), Yst(s, t));
   }
};

double integral(double f(double, double), const Triangle& tr, const QuadPoints& qp)
{
   const int order = qp.order;
   double result = 0;

   for(size_t i = 0; i < order; i++)
   {
      result += qp.w[i]* f(tr.Xst(qp.x[i], qp.y[i]), tr.Yst(qp.x[i], qp.y[i]));
   }

   return result * tr.GetArea();
}