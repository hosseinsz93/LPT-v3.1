#ifndef _GEOM_SPATIALINDEX_H
#define _GEOM_SPATIALINDEX_H


#include <memory>
#include <vector>

namespace Geom
{
  class BoundingBox;

  class Id;

  class SpatialIndexSource;

  /*! The SpatialIndex class is an abstraction used to represent
    a class of algorithms for storing objects in a spatial hierarchy
    such as an octree.  Any method that makes sense for any and all search
    trees (i.e. getOverlappingEntities()) should be defined as a pure virtual
    function in this base class so search trees can be swapped out when
    appropriate.

    The source argument passed into the SpatialIndex is cloned to deal with
    persistence issues.  In other words, the following code is valid:

    \code

    Rtree *createVertexTree(Geometry *geometry)
    {
      // create source object on local stack
      VertexSource source(geometry);
      // create new rtree object and return it
      return new Rtree(source);
      // source goes out of scope, but Rtree is valid as long as the 
      // Geometry object is valid
    }
    \endcode
  */
  class SpatialIndex
  {
  public:
    /*! Construct the spatial index.  The \a source argument passed 
      in must be a derived class of SpatialIndexSource and is used
      to look up coordinates, bounding boxes, and distances based on
      an Id.  Sources for standard elements (cells, faces, edges,
      vertices) are provided, but you can write your own source if
      needed also.

      A cloned copy of the source object is held by the SpatialIndex and
      deleted on its destruction.  Thus, it is safe to construct a SpatialIndex
      with a source defined in a local scope even if the index persists beyond
      that scope.
    */
    SpatialIndex(std::shared_ptr<SpatialIndexSource> const source);
    virtual ~SpatialIndex();
        
    /*! Retrieve the source object. */
    SpatialIndexSource const &getSource() const {return *_source;}

    /*! Get all entities who might be partially contained within the input
      bounding box.  The exact interpretation of this is up to the discretion
      of the derived class. */
    virtual void getOverlappingEntities(BoundingBox const &bbox,
                                        std::vector<Id> &entities) const = 0;


    /*! Get the bounding box containing all entities in the tree. */
    virtual BoundingBox getBoundingBox() const = 0;

  private:
    std::shared_ptr<SpatialIndexSource> const _source;
  };
}

#endif
