#ifndef _GEOM_SPHERE3
#define _GEOM_SPHERE3

#include "GeomBoundedShape3.h"

#include "GeomVector.h"

namespace Geom
{
  /*! A sphere in 3d. */
  class Sphere3 : public BoundedShape3
  {
  public:
    Sphere3(Vector3d const &center, double radius)
      : _center(center)
      , _radius(radius)
    {}

    Sphere3(Sphere3 const &rhs)
      : BoundedShape3(rhs),
	_center(rhs._center),
	_radius(rhs._radius)
    {}

    Sphere3 &operator=(Sphere3 const &rhs)
    {
      _center = rhs._center;
      _radius = rhs._radius;
      return *this;
    }

    Vector3d const &getCenter() const {
      return _center;
    }

    void setCenter(Vector3d const &center) {
      _center = center;
    }

    double const &getRadius() const {
      return _radius;
    }

    void setRadius(double radius) {
      _radius = radius;
    }

    Shape3 *clone() const {
      return new Sphere3(*this);
    }

    /*! Implement BoundedShape3::getBoundingBox() method. */
    BoundingBox getBoundingBox() const;

    /*! Implement BoundedShape3::intersectsBox() method. */
    bool intersectsBox(BoundingBox const &bbox) const;

    /*! Implement Shape::containsPoint() method. */
    bool containsPoint(Vector3d const &point) const {
      return (Vector3d(point - _center).mag2()
	      <= _radius * _radius);
    }

    bool isDegenerated() const;

  private:
    Vector3d _center;
    double _radius;
  };
}
    

#endif
