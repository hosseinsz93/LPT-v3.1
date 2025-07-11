#ifndef _GEOM_SPATIALINDEXDISTANCEQUERY_H
#define _GEOM_SPATIALINDEXDISTANCEQUERY_H

#include "GeomSpatialIndexQuery.h"

#include "GeomVector.h"

#include <limits> // for numeric_limits

namespace Geom
{
  class SpatialIndexSource;

  class DistanceSpatialIndexQuery : public SpatialIndexQuery
  {
  public:
    DistanceSpatialIndexQuery(SpatialIndexSource const &source,
			      Vector3d const &queryPoint,
			      double maxDistance = 
			      std::numeric_limits<double>::max(),
			      bool shrinkMaxDistance = false); 
    
    bool getDistance(Id const &id, double *distance) const;
    bool getDistance(Box3 const &bbox, double *distance) const;
    double convertDistance(double distance) const;

  private:
    SpatialIndexSource const &_source;
    Vector3d _queryPoint;
    mutable double _maxDistanceSq;
    bool _shrinkMaxDistance;
  };
}

#endif
