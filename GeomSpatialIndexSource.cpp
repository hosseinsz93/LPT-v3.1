#include "GeomSpatialIndexSource.h"

#include "GeomBoundingBox.h"
#include "GeomException.h"

#include "GeomId.h"

namespace Geom
{
  SpatialIndexSource::SpatialIndexSource()
  {

  }

  Vector3d SpatialIndexSource::getCoordinate(Id const &) const
  {
    throw Exception("getCoordinate not implemented for this "
                    "SpatialIndexSource");
    return Vector3d(0.);
  }

  BoundingBox SpatialIndexSource::getBoundingBox(Id const &)
    const
  {
    throw Exception("getBoundingBox not implemented for this "
                    "SpatialIndexSource");
    return BoundingBox();
  }

  double SpatialIndexSource::
  getDistanceSquared(Id const &id,
		     Vector3d const &p) const
  {
    throw Exception("getDistanceSquared not implemented for this "
                    "SpatialIndexSource");
  }

  // Added function so queries can be done on other geometric entities, e.g.
  // finding an edge closest to another edge.
  double SpatialIndexSource::
  getDistanceSquaredId(Id const & id,
                       MeshContainer3d const & queryContainer,
                       Id queryId) const
  {
    throw Exception("getDistanceSquaredId not implemented for this "
                    "SpatialIndexSource");
  }

  bool SpatialIndexSource::intersectsRay(Id const &/*id*/,
					 Ray3 const & /*ray*/,
					 double * /*distance*/) const 
  {
    throw Exception("intersectsRay not implemented for this "
                    "SpatialIndexSource");
    
  }
}
