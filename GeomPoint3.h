#ifndef _GEOM_POINT3
#define _GEOM_POINT3

#include "GeomBoundedShape3.h"

#include "GeomBoundingBox.h"
#include "GeomId.h"

namespace Geom
{
  class Point3 : public BoundedShape3
  {
  public:
    /** Construct an uninitialized point. */
    Point3()
    {}
    
    /** Initialize position from a Vector. */
    Point3(Vector3d const &position, Id const &id)
      : _position(position)
      , _id(id)
    {}

    /* note: use default copy constructor */

    /** Set/get the id */
    void setId(Id const &id){_id=id;}
    Id const &getId() const {return _id;}


    /** The position of the point. */
    Vector3d const &position() const { return _position; }
    /** The position of the point. */
    Vector3d       &position()       { return _position; }

    /** satisfy BoundedShape3 virtual interface. */
    BoundingBox getBoundingBox() const {
      return BoundingBox(_position);
    }

    bool intersectsBox(BoundingBox const &bbox) const {
      return bbox.containsPoint(_position);
    }

    /** satisfy Shape3 virtual interface. */
    Shape3 *clone() const {
      return new Point3(*this);
    }

    bool containsPoint(Vector3d const &point) const {
      return (Vector3d(point - position()).mag2() != 0);
    }

  private:
    Vector3d _position;
    Id _id;
  };

}

#endif
