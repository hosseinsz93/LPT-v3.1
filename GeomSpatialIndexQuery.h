#ifndef _GEOM_SPATIALINDEXQUERY_H
#define _GEOM_SPATIALINDEXQUERY_H

#include "GeomBox3.h"


namespace Geom
{
  class Id;

  /*! The SpatialIndexQuery is the base class for customized queries
    that can be passed into SpatialIndices like Rtrees.  Standard queries
    like distance from point and ray intersection are provided in the API
    of the Rtree, but customized queries can be defined here.  A custom query
    could be used, for instance, to query search tree for the closest entity
    to a given edge or face.

    All that needs to be defined for the query is the query distance for an
    object represented as an Id (this is the object actually stored in the 
    tree) and the query distance for a bounding box.  In order for searching
    to work properly, the query must be structured such that the query 
    distance for an object can never be less than the query distance for
    any box completely containing the object.

    The distance returned by the getDistance methods can actually be a 
    pseudo-distance for optimization reasons.  For instace, in a regular
    point distance query, it is typically cheaper to compute distance squared
    rather than distance.  In this case, the distance^2 can be returned by
    the getDistance functions, as long as the two getDistance functions return
    a consistent pseudo-distance.  In this case, the convertDistance method
    should also be overloaded to specity the conversion from pseudoDistance
    to real distance.
  */   
  class SpatialIndexQuery
  {
  public:
    SpatialIndexQuery() 
      : _hasQueryBox(false)
      {}

    virtual ~SpatialIndexQuery() {}

    /*! Get the query distance of an object stored in the search tree with
      Id \a id.  If the object is a valid candidate for the query (i.e.
      is of correct type, closer than maximum distance, etc.), the pseudo
      distance to the object is returned in argument \a distance and the 
      method should return true.  If the object is not valid, return false.
    */
    virtual bool getDistance(Id const &id,
			     double *distance) const = 0;

    /*! Get the query distance of a bounding box.  If the box is a valid 
      candidate for the query (i.e. closer than maximum distance), the 
      minimum pseudo distance to the box is returned in argument \a distance
      and the method should return true.  If the box is not valid, return 
      false.

      See notes in the class documentation about the relationship between
      this getDistance method and the one for objects.  For a typical point
      query, this means that the box getDistance method should return a distance
      of 0 if the query point is within the box and the closest distance 
      between the point and the box otherwise.
    */
    virtual bool getDistance(Box3 const &bbox,
			     double *distance) const = 0;

    /*! Convert between pseudo distance and real distance.  If the getDistance
      methods above both return something that is not the actual cartesian
      distance (i.e. distance squared), this method must be overloaded to
      specify conversion from pseudodistance to real distance.
    */
    virtual double convertDistance(double pseudoDistance) const {
      return pseudoDistance;
    }

    /*!
      Some types of queries can only be true if the object being tested
      intersects a fixed bounding box.  For example, point or ray queries
      with a finite maximum distance can only be satisfied for objects
      within that distance.  Since distance tests are generally slower
      than overlap tests in the rtree, queries like this can be faster
      if this overlap test is performed first.

      This method returns true if the query has a fixed box defined,
      false otherwise.
    */
    bool hasQueryBox() const {
      return _hasQueryBox;
    }

    /*!
     * See notes in hasQueryBox().  Return the query box for this query.
     */
    Box3 const &getQueryBox() const {
      return _queryBox;
    }

  protected:
    /*!
     * See notes in hasQueryBox().  Define the query box for this query
     */
    void setQueryBox(Box3 const &box) {
      _queryBox = box;
      _hasQueryBox = true;
    }

  private:
    bool _hasQueryBox;
    Box3 _queryBox;
  };
}

#endif
