#ifndef _SPATIALINDEXSURFS_H
#define _SPATIALINDEXSURFS_H

#include <memory>
#include <tuple>
#include <vector>

#include "GeomBoundingBox.h"
#include "GeomId.h"
#include "GeomRay3.h"
#include "GeomSpatialIndexSource.h"
#include "GeomVector.h"
#include "Surface.h"


struct UserCtx;

class SpatialIndexSurfs : public Geom::SpatialIndexSource
{
  public:
    SpatialIndexSurfs(std::shared_ptr< UserCtx > & _user,
                      std::shared_ptr< SpatialIndexSurfs > & _self);

    virtual ~SpatialIndexSurfs() {}

    /*! Clone the source to deal with persistence issues. All derived sources
      must implement this method to return a new, identical instance of the
      object.  See existing source objects for examples. */
    std::shared_ptr<Geom::SpatialIndexSource> clone() const;

    /*! Check if the given id corresponds to a point element (a single
      coordinate) or a region element (bounded by a finite bounding box).
    */
    bool isPoint(Geom::Id const &) const;

    /*! Get the coordinate value for the point.  This method will be
      called if isPoint() returns true for a given id. You do not have
      to implement this method in a derived class if it will never
      be called (i.e. isPoint() will always be false). */
    Vector3d getCoordinate(Geom::Id const &id) const;

    /*! Get the bounding box for the given object id.  This method will be
      called if isPoint() returns false for a given id. You do not have
      to implement this method in a derived class if it will never
      be called (i.e. isPoint() will always be true). */
    Geom::BoundingBox getBoundingBox(Geom::Id const &id) const;

    /*! Get the minimum squared distance between the test point \a p 
      and the given object id. 

      If either getDistanceSquared or getDistanceSquaredId is not implemented
      for a derived SpatialIndexSource, distance-based queries will not work.
    */
    double getDistanceSquared(Geom::Id const &id,
                              Vector3d const &p) const;

    /*! Get the minimum squared distance between the query entity defined by
      \a queryContainer and \a queryId and given object \a id. 

      If either getDistanceSquared or getDistanceSquaredId is not implemented
      for a derived SpatialIndexSource, distance-based queries will not work.
    */

    double getDistanceSquaredId(Geom::Id const & id,
                                Geom::Id queryId) const;

    /*! Determine if the object with given id intersects the specified ray.
      If so, return true and set /a distance to the distance along the ray 
      from its origin.  If not, return false.

      If this method is not implemented for a derived SpatialIndexSource,
      ray-based queries will not work.
    */
    bool intersectsRay(Geom::Id const &/*id*/,
                       Geom::Ray3 const & /*ray*/,
                       double * /*distance*/) const;
    
  private:
    std::shared_ptr< SpatialIndexSurfs > self;
    std::shared_ptr< UserCtx > user;
    std::shared_ptr< Surfaces > surfs;
};

#endif
