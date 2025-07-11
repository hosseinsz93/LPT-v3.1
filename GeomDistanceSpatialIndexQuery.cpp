#include "GeomDistanceSpatialIndexQuery.h"

#include "GeomSpatialIndexSource.h"

#include "GeomBox3.h"

#include <limits>

namespace Geom
{
  DistanceSpatialIndexQuery::
  DistanceSpatialIndexQuery(SpatialIndexSource const &source,
			    Vector3d const &queryPoint,
			    double maxDistance,
			    bool shrinkMaxDistance)
    :  _source(source),
       _queryPoint(queryPoint),
       _shrinkMaxDistance(shrinkMaxDistance)
  {
    // make sure not to overflow here
    if (maxDistance < sqrt(std::numeric_limits<double>::max()))
      {
	_maxDistanceSq = maxDistance*maxDistance;

	Geom::Vector<3, double> corner1(queryPoint(0) - maxDistance,
				        queryPoint(1) - maxDistance,
				        queryPoint(2) - maxDistance);
	Geom::Vector<3, double> corner2(queryPoint(0) + maxDistance,
				        queryPoint(1) + maxDistance,
				        queryPoint(2) + maxDistance);
	SpatialIndexQuery::setQueryBox(Box3(corner1, corner2));
      }
    else
      _maxDistanceSq = std::numeric_limits<double>::max();
  }

  bool DistanceSpatialIndexQuery::getDistance(Id const &id, 
					      double *distance) const
  {
    *distance = _source.getDistanceSquared(id, _queryPoint);
    if (*distance < _maxDistanceSq)
      {
	if (_shrinkMaxDistance)
	  _maxDistanceSq = *distance;
	return true;
      }
    return false;
  }

  bool DistanceSpatialIndexQuery::getDistance(Box3 const &bbox, 
					      double *distance) const
  {
    *distance = bbox.getDistanceSquared(_queryPoint);
    return (*distance < _maxDistanceSq);
  }
  
  double DistanceSpatialIndexQuery::convertDistance(double distance) const
  {
    return sqrt(distance);
  }
}

