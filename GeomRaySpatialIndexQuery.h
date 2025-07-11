#ifndef _GEOM_RAYSPATIALINDEXQUERY_H
#define _GEOM_RAYSPATIALINDEXQUERY_H

#include "GeomSpatialIndexQuery.h"

#include "GeomRay3.h"

#include <limits> // for numeric_limits

namespace Geom
{
  class SpatialIndexSource;

  class RaySpatialIndexQuery : public SpatialIndexQuery
  {
  public:
    RaySpatialIndexQuery(SpatialIndexSource const &source,
			 Ray3 * ray,
			 double maxDistance = 
			 std::numeric_limits<double>::max());
    bool getDistance(Id const &id, double *distance) const;
    bool getDistance(Box3 const &bbox, double *distance) const;
  private:
    SpatialIndexSource const &_source;
    Ray3 * _ray;
    double _maxDistance;
  };
}

#endif
