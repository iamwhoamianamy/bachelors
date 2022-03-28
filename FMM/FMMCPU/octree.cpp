#include "octree.hpp"
#include "math.hpp"
#include "spherical_harmonics.hpp"

Octree::Octree()
{
}

Octree::Octree(const Box& box, const size_t capacity) :
   _box(box), _capacity(capacity), _isSubdivided(false)
{

}

void Octree::insert(std::vector<Quadrature>& points)
{
   for(auto& point : points)
   {
      insert(point);
   }
}

void Octree::insert(Quadrature& point)
{
   if(!_box.contains(point.coordinates))
   {
      return;
   }

   if(_quadratures.size() < _capacity)
   {
      if(!_isSubdivided)
      {
         _quadratures.push_back(&point);
      }
      else
      {
         for(auto& child : _children)
         {
            child->insert(point);
         }
      }
   }
   else
   {
      subdivide();
      _isSubdivided = true;

      _quadratures.push_back(&point);

      for(auto p : _quadratures)
      {
         for(auto& child : _children)
         {
            child->insert(*p);
         }
      }

      _quadratures.clear();
   }
}

void Octree::subdivide()
{
   float x = _box.center.x;
   float y = _box.center.y;
   float z = _box.center.z;

   float w = _box.halfDimensions.x;
   float h = _box.halfDimensions.y;
   float d = _box.halfDimensions.z;

   Vector3 childrenHalfDimensions = { w / 2, h / 2, d / 2 };

   _children.reserve(8);

   _children.emplace_back(new Octree(Box({ x - w / 2, y - h / 2, z - d / 2 }, childrenHalfDimensions), _capacity));
   _children.emplace_back(new Octree(Box({ x + w / 2, y - h / 2, z - d / 2 }, childrenHalfDimensions), _capacity));
   _children.emplace_back(new Octree(Box({ x + w / 2, y + h / 2, z - d / 2 }, childrenHalfDimensions), _capacity));
   _children.emplace_back(new Octree(Box({ x - w / 2, y + h / 2, z - d / 2 }, childrenHalfDimensions), _capacity));

   _children.emplace_back(new Octree(Box({ x - w / 2, y - h / 2, z + d / 2 }, childrenHalfDimensions), _capacity));
   _children.emplace_back(new Octree(Box({ x + w / 2, y - h / 2, z + d / 2 }, childrenHalfDimensions), _capacity));
   _children.emplace_back(new Octree(Box({ x + w / 2, y + h / 2, z + d / 2 }, childrenHalfDimensions), _capacity));
   _children.emplace_back(new Octree(Box({ x - w / 2, y + h / 2, z + d / 2 }, childrenHalfDimensions), _capacity));
}

void Octree::quarry(const Box& range, std::vector<Quadrature*>& found)
{
   if(!this->_box.intersects(range))
   {
      return;
   }
   else
   {
      for(auto point : _quadratures)
      {
         if(range.contains(point->coordinates))
         {
            found.push_back(point);
         }
      }

      if(_isSubdivided)
      {
         for(auto& child : _children)
         {
            child->quarry(range, found);
         }
      }
   }
}

std::vector<Quadrature*> Octree::getAllQuadratures() const
{
   std::vector<Quadrature*> result;
   result.insert(result.end(), _quadratures.begin(), _quadratures.end());

   if(_isSubdivided)
   {
      for(auto child : _children)
      {
         auto childQudratures = child->getAllQuadratures();
         result.insert(result.end(), childQudratures.begin(), childQudratures.end());
      }
   }

   return result;
}

void Octree::calcLocalMultipolesWithoutTranslation(int n)
{
   _multipoleExpansion = math::calcIntegralContribution(getAllQuadratures(), n, _box.center);

   if(_isSubdivided)
   {
      for(auto& child : _children)
      {
         child->calcLocalMultipolesWithoutTranslation(n);
      }
   }
}

Vector3 Octree::calcA(const Vector3& point) const
{
   int n = _multipoleExpansion.data.size();

   Vector3 res;

   if(2 * _box.radius() < (point - _box.center).length())
   {
      auto irregularHarmonic = Harmonics::calcIrregularSolidHarmonics(n, point - _box.center);

      for(int l = 0; l < n; l++)
      {
         for(int m = -l; m <= l; m++)
         {
            res += _multipoleExpansion.getHarmonic(l, m) * irregularHarmonic.getHarmonic(l, m);
         }
      }
   }
   else
   {
      for(auto child : _children)
      {
         res += child->calcA(point);
      }
   }

   return res;
}

const Box& Octree::box() const
{
   return this->_box;
}

const bool Octree::isSubdivided() const
{
   return _isSubdivided;
}

const std::vector<Octree*> Octree::children() const
{
   return _children;
}

const HarmonicSeries<Vector3>& Octree::multipoleExpansion() const
{
   return _multipoleExpansion;
}

Octree::~Octree()
{
   if(_isSubdivided)
   {
      for(auto child : _children)
      {
         delete child;
      }
   }
}
