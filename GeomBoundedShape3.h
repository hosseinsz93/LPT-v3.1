#ifndef _GEOM_BOUNDEDSHAPE3
#define _GEOM_BOUNDEDSHAPE3

#include "GeomShape3.h"

namespace Geom
{
  class BoundingBox;

  /*! 
    The BoundedShape3 class is an abstract base class for 3-dimensional 
    geometrical shapes with finite bounds.  It provides abstract methods
    for retrieving the bounding box and testing for intersection with a
    bounding box (used for volume refinement sources).

    \todo This class currently deals with bounding boxes using a 
    BoundingBox object instead of a Box3 object -- these two objects
    could possibly be combined at some point...
  */
  class BoundedShape3 : public Shape3
  {
  public:
    virtual BoundingBox getBoundingBox() const  = 0;

    virtual bool intersectsBox(BoundingBox const &bbox) const = 0;

    virtual bool isDegenerate() const 
    { return false; }
  };
}

#endif
