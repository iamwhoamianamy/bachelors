#include "calculation_point_octree.hpp"
#include "math.hpp"
#include "integration.hpp"
#include "harmonics.hpp"
#include "multipole_translator.hpp"

#include <thrust/complex.h>

CalculationPointOctreeNode::CalculationPointOctreeNode() : 
   _parent(nullptr)
{
}

CalculationPointOctreeNode::CalculationPointOctreeNode(
   const Box& box,
   const size_t capacity,
   CalculationPointOctreeNode* parent) :
   _parent(parent), _box(box), _capacity(capacity)
{

}

void CalculationPointOctreeNode::insert(std::vector<Quadrature>& points)
{
   for(auto& point : points)
   {
      insert(point);
   }
}

void CalculationPointOctreeNode::insert(Quadrature& point)
{
   if(!_box.contains(point.coordinates))
   {
      return;
   }

   if(_quadratures.size() < _capacity)
   {
      if(!isSubdivided())
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

void CalculationPointOctreeNode::subdivide()
{
   float x = _box.center.x;
   float y = _box.center.y;
   float z = _box.center.z;

   float w = _box.halfDimensions.x;
   float h = _box.halfDimensions.y;
   float d = _box.halfDimensions.z;

   Vector3 childrenHalfDimensions = { w / 2, h / 2, d / 2 };

   _children.reserve(8);

   _children.emplace_back(new CalculationPointOctreeNode(Box({ x - w / 2, y - h / 2, z - d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new CalculationPointOctreeNode(Box({ x + w / 2, y - h / 2, z - d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new CalculationPointOctreeNode(Box({ x + w / 2, y + h / 2, z - d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new CalculationPointOctreeNode(Box({ x - w / 2, y + h / 2, z - d / 2 }, childrenHalfDimensions), _capacity, this));

   _children.emplace_back(new CalculationPointOctreeNode(Box({ x - w / 2, y - h / 2, z + d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new CalculationPointOctreeNode(Box({ x + w / 2, y - h / 2, z + d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new CalculationPointOctreeNode(Box({ x + w / 2, y + h / 2, z + d / 2 }, childrenHalfDimensions), _capacity, this));
   _children.emplace_back(new CalculationPointOctreeNode(Box({ x - w / 2, y + h / 2, z + d / 2 }, childrenHalfDimensions), _capacity, this));
}

void CalculationPointOctreeNode::quarry(const Box& range, std::vector<Quadrature*>& found)
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

      if(isSubdivided())
      {
         for(auto& child : _children)
         {
            child->quarry(range, found);
         }
      }
   }
}

std::vector<Quadrature*> CalculationPointOctreeNode::getAllQuadratures() const
{
   std::vector<Quadrature*> result;
   result.insert(result.end(), _quadratures.begin(), _quadratures.end());

   if(isSubdivided())
   {
      for(auto child : _children)
      {
         auto childQudratures = child->getAllQuadratures();
         result.insert(result.end(), childQudratures.begin(), childQudratures.end());
      }
   }

   return result;
}

void CalculationPointOctreeNode::calcLocalMultipolesWithoutTranslation(int n)
{
   if(!isLeafAndUseful())
   {
      _multipoleExpansion = math::calcIntegralContribution(getAllQuadratures(), n, _box.center);

      for(auto& child : _children)
      {
         child->calcLocalMultipolesWithoutTranslation(n);
      }
   }
}

void CalculationPointOctreeNode::calcLocalMultipolesWithComplexTranslation(int n)
{
   if(!isLeafAndUseful())
   {
      _multipoleExpansion = HarmonicSeries<Vector3>(n);

      for(auto child : _children)
      {
         child->calcLocalMultipolesWithComplexTranslation(n);
      }

      for(auto child : _children)
      {
         if(!child->_quadratures.empty() && !child->isSubdivided() ||
            child->_quadratures.empty() && child->isSubdivided())
            _multipoleExpansion.add(MultipoleTranslator::translateMultipoleWithComplex(
               child->multipoleExpansion(), child->box().center - _box.center));
      }
   }
}

void CalculationPointOctreeNode::calcLocalMultipolesWithRealTranslation(int n)
{
   if(!isLeafAndUseful())
   {
      _multipoleExpansion = HarmonicSeries<Vector3>(n);

      for(auto child : _children)
      {
         child->calcLocalMultipolesWithRealTranslation(n);
      }

      for(auto child : _children)
      {
         if(!child->_quadratures.empty() && !child->isSubdivided() ||
            child->_quadratures.empty() && child->isSubdivided())
            _multipoleExpansion.add(MultipoleTranslator::translateMultipoleWithReal(
               child->multipoleExpansion(), child->box().center - _box.center));
      }
   }
}

void CalculationPointOctreeNode::initAllMultipoleExpansions(size_t n)
{
   if(isSubdivided() || _quadratures.empty())
   {
      _multipoleExpansion = HarmonicSeries<Vector3>(n);

      for(auto child : _children)
      {
         child->initAllMultipoleExpansions(n);
      }
   }
}

void CalculationPointOctreeNode::calcLocalMultipolesAtLeaves(size_t n)
{
   if(isLeafAndUseful())
   {
      _multipoleExpansion = math::calcIntegralContribution(
         _quadratures, n, _box.center);
   }
   else
   {
      for(auto child : _children)
      {
         child->calcLocalMultipolesAtLeaves(n);
      }
   }
}

Vector3 CalculationPointOctreeNode::calcA(const Vector3& point) const
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

Vector3 CalculationPointOctreeNode::caclRot(const Vector3& point) const
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

size_t CalculationPointOctreeNode::getAllNodeCount() const
{
   size_t count = 1;

   if(isSubdivided())
   {
      count += 8;

      for(auto child : _children)
      {
         count += child->getAllNodeCount();
      }
   }
   
   return count;
}

const Box& CalculationPointOctreeNode::box() const
{
   return this->_box;
}

bool CalculationPointOctreeNode::isSubdivided() const
{
   return _children.size();
}

CalculationPointOctreeNode* CalculationPointOctreeNode::parent()
{
   return _parent;
}

CalculationPointOctreeNode* CalculationPointOctreeNode::parent() const
{
   return _parent;
}

std::vector<CalculationPointOctreeNode*>& CalculationPointOctreeNode::children()
{
   return _children;
}

const std::vector<CalculationPointOctreeNode*>& CalculationPointOctreeNode::children() const
{
   return _children;
}

std::vector<Quadrature*> CalculationPointOctreeNode::quadratures() const
{
   return _quadratures;
}

HarmonicSeries<Vector3>& CalculationPointOctreeNode::multipoleExpansion()
{
   return _multipoleExpansion;
}

const HarmonicSeries<Vector3>& CalculationPointOctreeNode::multipoleExpansion() const
{
   return _multipoleExpansion;
}

CalculationPointOctreeNode::~CalculationPointOctreeNode()
{
   if(isSubdivided())
   {
      for(auto child : _children)
      {
         delete child;
      }
   }
}

bool CalculationPointOctreeNode::isLeafAndUseful() const
{
   return !isSubdivided() && !_quadratures.empty();
}
