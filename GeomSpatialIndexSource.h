#ifndef _GEOM_SPATIALINDEXSOURCE_H
#define _GEOM_SPATIALINDEXSOURCE_H

#include <memory>

#include "GeomVector.h"

namespace Geom
{
  class BoundingBox;

  class Ray3;
  class Point3;

  class Id;
  class MeshContainer3d;

  /*! The SpatialIndexSource is an abstraction between a spatial index
    (i.e. an octree) and its data.  It provides functions to look up
    coordinate, bounding box, and distance from a query point. 

    The SpatialIndexSource uses a clone idiom to provide a pass-by-value
    feel while remaining polymorphic.  Thus, derived sources should be 
    lightweight objects, only holding on to pointers to other data (i.e. 
    TopologyData, etc.).
  */
  class SpatialIndexSource
  {
  public:
    SpatialIndexSource();
    virtual ~SpatialIndexSource() {}

    /*! Clone the source to deal with persistence issues. All derived sources
      must implement this method to return a new, identical instance of the
      object.  See existing source objects for examples. */
      virtual std::shared_ptr<SpatialIndexSource> clone() const = 0;

    /*! Check if the given id corresponds to a point element (a single
      coordinate) or a region element (bounded by a finite bounding box).
    */
    virtual bool isPoint(Id const &) const  = 0;

    /*! Get the coordinate value for the point.  This method will be
      called if isPoint() returns true for a given id. You do not have
      to implement this method in a derived class if it will never
      be called (i.e. isPoint() will always be false). */
    virtual Vector3d getCoordinate(Id const &id) const;

    /*! Get the bounding box for the given object id.  This method will be
      called if isPoint() returns false for a given id. You do not have
      to implement this method in a derived class if it will never
      be called (i.e. isPoint() will always be true). */
    virtual BoundingBox getBoundingBox(Id const &id) const;

    /*! Get the minimum squared distance between the test point \a p 
      and the given object id. 

      If either getDistanceSquared or getDistanceSquaredId is not implemented
      for a derived SpatialIndexSource, distance-based queries will not work.
    */
    virtual double getDistanceSquared(Id const &id,
                                      Vector3d const &p) const;

    /*! Get the minimum squared distance between the query entity defined by
      \a queryContainer and \a queryId and given object \a id. 

      If either getDistanceSquared or getDistanceSquaredId is not implemented
      for a derived SpatialIndexSource, distance-based queries will not work.
    */

    virtual double
    getDistanceSquaredId(Id const & id,
                         MeshContainer3d const & queryContainer,
                         Id queryId) const;

    /*! Determine if the object with given id intersects the specified ray.
      If so, return true and set /a distance to the distance along the ray 
      from its origin.  If not, return false.

      If this method is not implemented for a derived SpatialIndexSource,
      ray-based queries will not work.
    */
    virtual bool intersectsRay(Id const &/*id*/,
                               Ray3 const & /*ray*/,
                               double * /*distance*/) const;
  };
}

#endif  
