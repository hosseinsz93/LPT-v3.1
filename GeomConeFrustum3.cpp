#include "GeomConeFrustum3.h"

#include "GeomBox3.h"
#include "GeomPoint3.h"

#include "GeomUtility.h"


namespace Geom
{
  ConeFrustum3::ConeFrustum3(Vector3d const &center1, double radius1,
			     Vector3d const &center2, double radius2)
    : _centerLineLength(0)
    , _hasBoundingBox(false)
    , _hasConeVertex(false)
  {
    if (radius1 < 0)
      radius1 *= -1;
    if (radius2 < 0)
      radius2 *= -1;
    
    if (radius1 > radius2)
      {
	_radius1 = radius2;
	_center1 = center2;
	_radius2 = radius1;
	_center2 = center1;
      }
    else 
      {
	_radius1 = radius1;
	_center1 = center1;
	_radius2 = radius2;
	_center2 = center2;
      }
    _centerLineLength = _center1.distance(_center2);
  }

  ConeFrustum3::ConeFrustum3(ConeFrustum3 const &rhs)
    : BoundedShape3(rhs)
    , _center1(rhs._center1)
    , _radius1(rhs._radius1)
    , _center2(rhs._center2)
    , _radius2(rhs._radius2)
    , _centerLineLength(rhs._centerLineLength)
    , _hasBoundingBox(rhs._hasBoundingBox)
    , _boundingBox(rhs._boundingBox)
    , _hasConeVertex(rhs._hasConeVertex)
    , _coneVertex(rhs._coneVertex)
  {}

  ConeFrustum3 & ConeFrustum3::operator=(ConeFrustum3 const &rhs)
  {
    if (&rhs == this)
      return *this;

    _center1 = rhs._center1; 
    _radius1 = rhs._radius1;
    _center2 = rhs._center2; 
    _radius2 = rhs._radius2;
    _hasBoundingBox = rhs._hasBoundingBox;
    _boundingBox = rhs._boundingBox;
    _hasConeVertex = rhs._hasConeVertex;
    _coneVertex = rhs._coneVertex;
    return *this;
  }

  bool ConeFrustum3::intersectsBox(BoundingBox const &bbox) const
  {
    return Box3(bbox.min(), bbox.max()).testIntersection(*this);
  }

  bool ConeFrustum3::containsPoint(Vector3d const &point) const
  {
    return testIntersection(point);
  }
  
  void ConeFrustum3::computeBoundingBox() const
  {
    Vector3d extent(_center2-_center1);
    double Rextent = extent.length();
    Vector3d normalExtent = extent;
    normalExtent.normalize();
    Vector3d khat(0,0,1);
    Vector3d jhat(0,1,0);
    Vector3d ihat(1,0,0);
    double cosz = normalExtent.dot(khat);
    double cosy = normalExtent.dot(jhat);
    double cosx = normalExtent.dot(ihat);
    
    double sinz = sqrt(1.0-cosz*cosz);
    double siny = sqrt(1.0-cosy*cosy);
    double sinx = sqrt(1.0-cosx*cosx);
    
    double z1,z2,z3,y1,y2,y3,x1,x2,x3;
    double minx,miny,minz,maxx,maxy,maxz;
    
    z1 = Rextent*cosz + _radius2*sinz;
    z2 = Rextent*cosz - _radius2*sinz;
    y1 = Rextent*cosy + _radius2*siny;
    y2 = Rextent*cosy - _radius2*siny;
    x1 = Rextent*cosx + _radius2*sinx;
    x2 = Rextent*cosx - _radius2*sinx;
    
    z3 = _radius1*sinz;
    y3 = _radius1*siny;
    x3 = _radius1*sinx;
    if (z3 < -z3)
      {
	minz = z3;
	maxz = -z3;
      }
    else 
      {
	minz = -z3;
	maxz = z3;
      }
    if (x3 < -x3)
      {
	minx = x3;
	maxx = -x3;
      }
    else
      {
	minx = -x3;
	maxx = x3;
      }
    if (y3 < -y3)
      {
	miny = y3;
	maxy = -y3;
      }
    else
      {
	miny = -y3;
	maxy = y3;
      }
    minz = std::min(std::min(z1,z2),minz);
    miny = std::min(std::min(y1,y2),miny);
    minx = std::min(std::min(x1,x2),minx);
    maxz = std::max(std::max(z1,z2),maxz);
    maxy = std::max(std::max(y1,y2),maxy);
    maxx = std::max(std::max(x1,x2),maxx);

    _boundingBox = 
      BoundingBox(Vector3d(minx,miny,minz)+_center1,
                  Vector3d(maxx,maxy,maxz)+_center1);
    _hasBoundingBox = true;
  }

  double ConeFrustum3::getConeAngle() const
  {
    return atan2((_radius2-_radius1),
                  _center2.distance(_center1));
  }

  bool ConeFrustum3::testIntersection(Vector3d const &coord) const
  {
    // precompute some vectors
    Vector3d axis = this->getAxisVector();
    double axisMag = axis.normalize();
    if (axisMag <= ALMOST_ZERO)
      return false;

    Vector3d pointToConeStart(coord - this->getCenter1());

    // compute distance along axis between point and frustum
    double distAlongAxis = axis.dot(pointToConeStart);
    if (distAlongAxis < 0.0 || distAlongAxis > axisMag)
      return false; // point beyond axis of frustum

    // compute distance from point to frustum axis
    double distanceFromAxisSq = 
      pointToConeStart.mag2() - distAlongAxis * distAlongAxis;

    // compute frustum radius at point on axis
    double radiusAtPoint = this->getRadius1() + (this->getRadius2() - this->getRadius1()) * distAlongAxis / axisMag;

    return (distanceFromAxisSq <= radiusAtPoint * radiusAtPoint);
  }

  bool ConeFrustum3::isDegenerate() const
  {
    if ( (getRadius1() < ALMOST_ZERO) &&  (getRadius2() < ALMOST_ZERO))
      return true;

    if ( _centerLineLength < ALMOST_ZERO)
      return true;

    //is NOT degenerated cylinder
    return false;
  }
}   
