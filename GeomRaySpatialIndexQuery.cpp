#include "GeomRaySpatialIndexQuery.h"

#include "GeomSpatialIndexSource.h"

#include "GeomBox3.h"

namespace Geom
{  
  RaySpatialIndexQuery::RaySpatialIndexQuery(SpatialIndexSource const &source,
					     Ray3 * ray,
					     double maxDistance) :
    _source(source),
    _ray(ray),
    _maxDistance(maxDistance)
  {
    if (_maxDistance < std::numeric_limits<double>::max())
      {
	// maximum distance has been specified.  Compute a bounding
	// box containing the ray origin and it's "endpoint" based
	// on the distance and ray direction and store this
	// as the query box of this query to speed up searching.
	Vector3d const &origin = ray->getOrigin();
	Vector3d othercorner(origin);
	othercorner += ray->getDirection() * _maxDistance;

	SpatialIndexQuery::setQueryBox(Box3(origin, othercorner));
      }
  }

  bool RaySpatialIndexQuery::getDistance(Id const &id, 
					 double *distance) const
  {
    // rewrite this once source class is cleaned up
    if (!_source.intersectsRay(id, *_ray, distance))
      return false;
    return (*distance < _maxDistance);
  }
  
  bool RaySpatialIndexQuery::getDistance(Box3 const &box, 
					 double *distance) const
  {
    if (!box.testIntersection(*_ray, distance))
      return false;
    return (*distance < _maxDistance);
  }
  
}
