#include "GeomBox3.h"
#include "GeomSphere3.h"

#include "GeomBoundingBox.h"
#include "GeomUtility.h"

namespace Geom
{
  BoundingBox Sphere3::getBoundingBox() const
  {
    return BoundingBox(Vector3d(_center - _radius),
                       Vector3d(_center + _radius));
  }

  bool Sphere3::intersectsBox(BoundingBox const &bbox) const
  {
    return Box3(bbox.min(), bbox.max()).testIntersection(*this);
  }

  bool Sphere3::isDegenerated() const
  {
    if (getRadius() <= ALMOST_ZERO)
      return true;
    return false;
  }
}
