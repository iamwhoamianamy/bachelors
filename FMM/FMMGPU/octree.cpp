#include "octree.hpp"
#include "math.hpp"
#include "integration.hpp"
#include "harmonics.hpp"
#include <complex>

OctreeNode::OctreeNode() : 
   _parent(nullptr)
{
}

OctreeNode::OctreeNode(const Box& box, const size_t capacity, OctreeNode* parent) :
   _parent(parent), _box(box), _capacity(capacity), _isSubdivided(false)
{

}

void OctreeNode::insert(std::vector<Quadrature>& points)
{
   for(auto& point : points)
   {
      insert(point);
   }
}

void OctreeNode::insert(Quadrature& point)
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

void OctreeNode::subdivide()
{
   float x = _box.center.x;
   float y = _box.center.y;
   float z = _box.center.z;

   float w = _box.halfDimensions.x;
   float h = _box.halfDimensions.y;
   float d = _box.halfDimensions.z;

   Vector3 childrenHalfDimensions = { w / 2, h / 2, d / 2 };

   _children.reserve(8);

   _children.emplace_back(new OctreeNode(Box({ x - w / 2, y - h / 2, z - d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new OctreeNode(Box({ x + w / 2, y - h / 2, z - d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new OctreeNode(Box({ x + w / 2, y + h / 2, z - d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new OctreeNode(Box({ x - w / 2, y + h / 2, z - d / 2 }, childrenHalfDimensions), _capacity, this));

   _children.emplace_back(new OctreeNode(Box({ x - w / 2, y - h / 2, z + d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new OctreeNode(Box({ x + w / 2, y - h / 2, z + d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new OctreeNode(Box({ x + w / 2, y + h / 2, z + d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new OctreeNode(Box({ x - w / 2, y + h / 2, z + d / 2 }, childrenHalfDimensions), _capacity, this));
}

void OctreeNode::quarry(const Box& range, std::vector<Quadrature*>& found)
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

std::vector<Quadrature*> OctreeNode::getAllQuadratures() const
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

void OctreeNode::calcLocalMultipolesWithoutTranslation(int n)
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

void OctreeNode::calcLocalMultipolesWithComplexTranslation(int n)
{
   if(!_isSubdivided && !_quadratures.empty())
   {
      _multipoleExpansion = math::calcIntegralContribution(_quadratures, n, _box.center);
   }
   else
   {
      _multipoleExpansion = HarmonicSeries<Vector3>(n);

      for(auto child : _children)
      {
         child->calcLocalMultipolesWithComplexTranslation(n);
      }

      for(auto child : _children)
      {
         if(!child->_quadratures.empty() && !child->_isSubdivided ||
            child->_quadratures.empty() && child->_isSubdivided)
            _multipoleExpansion.add(Harmonics::translateWithComplex(
               child->multipoleExpansion(), child->box().center - _box.center));
      }
   }
}

void OctreeNode::calcLocalMultipolesWithRealTranslation(int n)
{
   if(!_isSubdivided && !_quadratures.empty())
   {
      _multipoleExpansion = math::calcIntegralContribution(_quadratures, n, _box.center);
   }
   else
   {
      _multipoleExpansion = HarmonicSeries<Vector3>(n);

      for(auto child : _children)
      {
         child->calcLocalMultipolesWithRealTranslation(n);
      }

      for(auto child : _children)
      {
         if(!child->_quadratures.empty() && !child->_isSubdivided ||
             child->_quadratures.empty() && child->_isSubdivided)
            _multipoleExpansion.add(Harmonics::translateWithReal(
               child->multipoleExpansion(), child->box().center - _box.center));
      }
   }
}

void OctreeNode::initAllMultipoleExpansions(size_t n)
{
   if(_isSubdivided || _quadratures.empty())
   {
      _multipoleExpansion = HarmonicSeries<Vector3>(n);

      for(auto child : _children)
      {
         child->initAllMultipoleExpansions(n);
      }
   }
}

Vector3 OctreeNode::calcA(const Vector3& point) const
{
   int n = _multipoleExpansion.order();

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

Vector3 OctreeNode::caclRot(const Vector3& point) const
{
   int n = _multipoleExpansion.order();
   real eps = 1e-3;

   Vector3 res;

   if(2 * _box.radius() < (point - _box.center).length())
   {
      auto hx1 = Harmonics::calcIrregularSolidHarmonics(n, point - _box.center + 
                                                        Vector3::xAxis() * eps);
      auto hx2 = Harmonics::calcIrregularSolidHarmonics(n, point - _box.center -
                                                        Vector3::xAxis() * eps);

      auto hy1 = Harmonics::calcIrregularSolidHarmonics(n, point - _box.center +
                                                        Vector3::yAxis() * eps);
      auto hy2 = Harmonics::calcIrregularSolidHarmonics(n, point - _box.center -
                                                        Vector3::yAxis() * eps);

      auto hz1 = Harmonics::calcIrregularSolidHarmonics(n, point - _box.center +
                                                        Vector3::zAxis() * eps);
      auto hz2 = Harmonics::calcIrregularSolidHarmonics(n, point - _box.center -
                                                        Vector3::zAxis() * eps);

      hx1.subtract(hx2);
      hy1.subtract(hy2);
      hz1.subtract(hz2);

      Vector3 tempRes;

      for(int l = 0; l < n; l++)
      {
         for(int m = -l; m <= l; m++)
         {
            tempRes.x += _multipoleExpansion.getHarmonic(l, m).z * hy1.getHarmonic(l, m) -
                         _multipoleExpansion.getHarmonic(l, m).y * hz1.getHarmonic(l, m);

            tempRes.y += _multipoleExpansion.getHarmonic(l, m).x * hz1.getHarmonic(l, m) -
                         _multipoleExpansion.getHarmonic(l, m).z * hx1.getHarmonic(l, m);

            tempRes.z += _multipoleExpansion.getHarmonic(l, m).y * hx1.getHarmonic(l, m) -
                         _multipoleExpansion.getHarmonic(l, m).x * hy1.getHarmonic(l, m);
         }
      }

      res += tempRes / (2 * eps);
   }
   else
   {
      for(auto child : _children)
      {
         res += child->caclRot(point);
      }
   }

   return res;
}

const Box& OctreeNode::box() const
{
   return this->_box;
}

const bool OctreeNode::isSubdivided() const
{
   return _isSubdivided;
}

OctreeNode* OctreeNode::parent()
{
   return _parent;
}

const OctreeNode* OctreeNode::parent() const
{
   return _parent;
}

const std::vector<OctreeNode*> OctreeNode::children() const
{
   return _children;
}

const std::vector<Quadrature*> OctreeNode::quadratures() const
{
   return _quadratures;
}

HarmonicSeries<Vector3>& OctreeNode::multipoleExpansion()
{
   return _multipoleExpansion;
}

const HarmonicSeries<Vector3>& OctreeNode::multipoleExpansion() const
{
   return _multipoleExpansion;
}

OctreeNode::~OctreeNode()
{
   if(_isSubdivided)
   {
      for(auto child : _children)
      {
         delete child;
      }
   }
}
