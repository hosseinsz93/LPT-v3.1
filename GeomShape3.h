#ifndef _GEOM_SHAPE3_H
#define _GEOM_SHAPE3_H

#include "GeomVector.h"

namespace Geom
{
  /*!
    The Shape3 class is an abstract base class for 3-dimensional geometric 
    shapes.  All 3d shapes (Box3, Sphere3, Line3, etc.) will derive from 
    this class.

    Any algorithms that need to operate on geometric shapes but do so
    without depending on concrete instances of the shapes can refer to the
    base class.  Feel free to add virtual methods to this class as needed
    to support such algorithms.
  */
  class Shape3
  {
  public:
    virtual ~Shape3() {}

    /*! Clone method to avoid slicing problems. */
    virtual Shape3 *clone() const = 0;

    /*! Return true if the shape contains the given point. */
    virtual bool containsPoint(Vector3d const &point) const = 0;
  };
}

#endif
