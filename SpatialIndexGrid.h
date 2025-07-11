#ifndef _SPATIALINDEXGRID_H
#define _SPATIALINDEXGRID_H

#include <memory>
#include <vector>

#include "GeomBoundingBox.h"
#include "GeomId.h"
#include "GeomRay3.h"
#include "GeomSpatialIndexSource.h"
#include "GeomVector.h"

#include "variables.h"


struct CmpntsStruct;

class SpatialIndexGrid : public Geom::SpatialIndexSource
{
  public:
    enum SidednessEnum
    {
      INTERNAL,
      ONBOUNDARY
    };

    struct Sidedness
    {
      Sidedness(SidednessEnum s) : sidedness(s) {}
      bool operator==(SidednessEnum const & s) const
      {
        return sidedness == s;
      }
      SidednessEnum sidedness;
    };
    
    SpatialIndexGrid(std::shared_ptr< UserCtx > & _user,
                     std::shared_ptr< SpatialIndexGrid > & _self);

    virtual ~SpatialIndexGrid() {}

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


    
    static void getCellCoords(int const i, int const j, int const k,
                              Cmpnts const * const * const * const coor,
                              std::vector<Vector3d> & fcoor);

    static void getVertCoords(int vertId, int const i, int const j, int const k,
                              Cmpnts const * const * const * const coor,
                              Vector3d & fcoor);

    static void getFaceCoords(int faceId, int const i, int const j, int const k,
                              Cmpnts const * const * const * const coor,
                              std::vector<Vector3d> & fcoor);
    
    void indexRanges(PetscInt & _xs, PetscInt & _xe,
                     PetscInt & _ys, PetscInt & _ye,
                     PetscInt & _zs, PetscInt & _ze)
    {
      _xs = xs;
      _xe = xe;
      _ys = ys;
      _ye = ye;
      _zs = zs;
      _ze = ze;
    }

    static Sidedness getCellSidedness(int i, int j, int k, int & iface);

    double getMaxDelta() const
    {
      return maxDelta;
    }

    typedef const int (&fidxs)[3];
    static fidxs getFaceIds(int fidx)
    {
      return n[fidx];
    }
    
    
  private:
    std::shared_ptr< SpatialIndexGrid > self;
    std::shared_ptr< UserCtx > user;

    static int ni, nj, nk;

    static PetscInt xs, xe;
    static PetscInt ys, ye;
    static PetscInt zs, ze;
    static PetscInt mx, my, mz;

    double maxDelta;

    // When running the intersectsRay, it knows which face is the one
    // pierced for the return value.  This information is useful since
    // the face corresponds to the grid i,j,k place that was pierced and
    // tells which plane the particle pierced to exit the grid.

    Cmpnts const * const * const * const coor;

    // Face definitions for a grid cell. Node offsets
    static constexpr int f[6][4]
    { { 1, 2, 6, 5 },    // 0
      { 0, 4, 7, 3 },    // 1
      { 3, 7, 6, 2 },    // 2
      { 0, 1, 5, 4 },    // 3
      { 4, 5, 6, 7 },    // 4
      { 0, 3, 2, 1 } };  // 5

    
    // Node definitions for a grid cell. Offsets from i,j,k.
    static constexpr int n[8][3]
    { { 0, 0, 0 },    // 0
      { 1, 0, 0 },    // 1
      { 1, 1, 0 },    // 2
      { 0, 1, 0 },    // 3
      { 0, 0, 1 },    // 4
      { 1, 0, 1 },    // 5
      { 1, 1, 1 },    // 6
      { 0, 1, 1 } };  // 7

    // Edge definition for a grid cell. Offsets from i,j,k.
    static constexpr int e[12][2][3]
    { { { 0, 0, 0 }, { 1, 0, 0 } },    // 0-1
      { { 1, 0, 0 }, { 1, 1, 0 } },    // 1-2
      { { 1, 1, 0 }, { 0, 1, 0 } },    // 2-3
      { { 0, 1, 0 }, { 0, 0, 0 } },    // 3-0

      { { 0, 0, 0 }, { 0, 0, 1 } },    // 0-4
      { { 1, 0, 0 }, { 1, 0, 1 } },    // 1-5
      { { 1, 1, 0 }, { 1, 1, 1 } },    // 2-6
      { { 0, 1, 0 }, { 0, 1, 1 } },    // 3-7
        
      { { 0, 0, 1 }, { 1, 0, 1 } },    // 4-5
      { { 1, 0, 1 }, { 1, 1, 1 } },    // 5-6
      { { 1, 1, 1 }, { 0, 1, 1 } },    // 6-7
      { { 0, 1, 1 }, { 0, 0, 1 } } };  // 7-4
};

#endif
