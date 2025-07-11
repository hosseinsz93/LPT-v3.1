#ifndef _GEOM_RTREE_H
#define _GEOM_RTREE_H

#include "GeomSpatialIndex.h"
#include "GeomVector.h"

#include <limits>
#include <memory>
#include <vector>

// #define DEBUG_RTREE_WRO

template <int N, typename T> class Vector;

namespace Geom
{
  class Subset;

  class Point3;
  class Ray3;

  class SpatialIndexQuery;


  /*! The Rtree class is an implementation of an R-tree spatial index tree.
    The R-tree was originally proposed in the paper "R-trees: A Dynamic
    Index Structure for Spatial Searching" by Guttman.  

    Basically, the R-tree is a search tree with the following properties:
    * Height-balanced, i.e. all leaf nodes are at the same depth.
    * All nodes (leaf and non-leaf) contain a fixed maximum number 
    (currently hardcoded to N=8) of children.  This allows the
    leaf depth of the tree to be relatively low (d ~= log_N(n) where n is 
    the number of leaf entries) which uses less total nodes for a given
    number of leaf entries. 
    * All nodes store the Minimum Bounding Rectangle (MBR) of all their 
    children's MBRs, recursively built up from the MBR of the actual
    leaf entries.
    * Overlapping, i.e. All leaf entries are fully contained in a single 
    parent but data in nodes can overlap

    A modified version of a published "bulk-loading" algorithm is implemented
    in the Rtree.  For details, see the notes in the .cpp file.  This algorithm
    is invoked in the addSubset method and is used to quickly and efficiently
    pack the R-tree given a list of leaf entities. 

    \todo Add incremental adding/removal of elements
    \todo Add addition of subsets (following "Small Tree/Large Tree" papers)
  */
  class Rtree : public SpatialIndex
  {
  public:
    /*! Construct the R-tree with given input source. */
    Rtree(std::shared_ptr<SpatialIndexSource> const & source);
    ~Rtree();

    /*! Query the R-tree for entries overlapping the input 
      bounding box.  Any entry whos bounding box overlaps the input
      box will be pushed onto the back of the \a entities vector.  */
    void getOverlappingEntities(BoundingBox const &bbox,
                                std::vector<Id> &entities) const;
    /*! Query the R-tree for entries overlapping the input 
      bounding box.  Entries are not returned just a bool to say if there
      are at least one  */
    bool checkOverlappingEntities(BoundingBox const &bbox) const;
    
    /*! Add a subset of elements to the R-tree.  APPENDING OF ADDITIONAL
      SUBSETS TO A NON-EMPTY R-TREE NOT YET IMPLEMENTED! */
    void addSubset(Subset const &subset);

    /*! Implement this method from base SpatialIndex. */
    BoundingBox getBoundingBox() const;

    /*! Clear the r-tree, destroying all nodes. */
    void clear();

    /*! Begin a point query.  Subsequent calls to query() will return 
      objects in increasing order of distance from the query point.  If
      maxDistance is set to a valid value, the query will terminate after
      no possible candidates closer to this distance can be found.
    */
    void beginPointQuery(Vector3d const &point,
			 double maxDistance = 
			 std::numeric_limits<double>::max()) const;

    /*! Begin a ray query.  Subsequent calls to query() will return
      objects that are pierced by the ray in increasing order of distance
      from its origin. 
      If maxDistance is set to a valid value, the query will terminate after
      no possible candidates closer to this distance can be found.
    */
    void beginRayQuery(Ray3 const &ray,
		       double maxDistance =  
		       std::numeric_limits<double>::max()) const;

    /*! Begin a customized query.  See SpatialIndexQuery for details on 
      constructing your own query functions.

      Ownership of the query pointer is transfered to the Rtree -- it will
      be deleted upon a call to endQuery().
    */
    void beginQuery(SpatialIndexQuery *query) const;
    
    /*! Process the current query function.  Must be called after one of the
      begin*Query functions and before a call to endQuery.  Returns the Id of
      the next object as defined by the query function or an invalid Id if
      no more objects can be found.  The actual distance (as defined by the
      query function) is placed in the return argument object_distance if
      non-null. */
    Id query(double *object_distance = 0) const;

    /*! End the current query.  Must be called after begin*Query. */
    void endQuery() const;

    /*! An optimized query that returns the single closest entity to
      the given point.  This is a bit faster than using the distance browsing 
      approach, but obviously provides less flexibility.

      \a point The query point
      \a object_distance If non-null, the actual distance from the input point 
      to the found object will be placed in here.
      \a maxDistance If specified, objects farther than this distance won't be
      considered during the search.  In some cases, this will speed up the
      search.

      \return The Id of the closest entity to the point, if one exists within
      specified maxDistance.  Otherwise, an invalid Id will be returned.
    */
    Id getClosestEntity(Vector3d const &point,
                        double *object_distance = 0,
                        double maxDistance = 
                        std::numeric_limits<double>::max()) const;    

    /*!
      Recompute all cached bounding boxes in the rtree without rebuilding
      its topology.  This keeps the hierarchy of the tree the same but will
      force recomputation of all leaf node bounding boxes based on the current
      position of elements in the input source.  These bounding box changes
      will then be propagated all the way up to the root.

      This method should be considerably faster to call than clearing out the
      entire rtree and rebuilding it from scratch.  As the motion of elements
      becomes large, though, the quality of the tree (and therefore the 
      efficiency of searches) will tend to decrease.
    */
    void recomputeAllBoundingBoxes();

    /*!
      Recompute the bounding box of the leaf node containing the given Id
      and propagate this modified bounding box up to the root.  See notes
      and warnings on the previous method for more details.
      WARNING: This is O(n).  Calling this in a loop will be very slow.
      Better is to collect all faces in a Subset and do all of them at once.
    */
    void recomputeBoundingBoxesContaining(Id const &id);

    /*!
      Recompute the bounding box of all leaf nodes containing elements in
      the given subset and propagate this modified bounding box up to the root.
      See notes and warnings on the previous methods for more details.
      Note: This is O(n) so time scales with size of Rtree, not size of Subset.
    */
    void recomputeBoundingBoxesContaining(Subset const &subset);

    /*! 
      add an element to the rtree.  This method should not be used
      to create the rtree from a large set of elements "addSubset"
      will do a better job of packing the tree and make a more efficient
      tree.  This is used when faces are added during a process.  If 
      a huge number of elements are added you may decide that the 
      time for rebuilding the tree is worth the performance savings of
      a more efficient tree.  
    */
    bool addElement(Id item);

    /*!
      remove an element from the rtree.  NOTE: This method requires that the
      element is still defined in the source and that it's bounding box is
      the same as that cached in the rtree.  In other words, you should never
      call this method after deleting the element from the data feeding the 
      source or after changing the bounding box of the object.  

      If you need to do this, call removeInvalidElement below.
    */
    bool removeElement(Id item);

    /*!
      Remove an element from the rtree by searching for it purely 
      topologically, i.e. w/o using the bounding box.  Use this when the
      constraints of the removeElement call cannot be satisfied.
    */
    bool removeInvalidElement(Id item);

    /*
      Remove a subset of elements from the tree.  This method will traverse
      the leaves of the tree looking for elements in the input subset to 
      remove.    No assumptions on the validity of the elements are made,
      i.e. they could have already been removed from the source ala 
      removeInvalidElement.

      All elements in the subset are removed from the tree.   The current
      implementation follows the simple but extremely fast "free-at-empty"
      technique described in "Lazy Deletion Methods in R-Trees" by A. 
      Nanopoulos et.al.   This will ensure that any empty nodes (leaf or 
      otherwise) will be pruned from the tree while underfull but not-empty
      nodes will remain.    A subsequent step of pruning underfull nodes could
      also be implemented to improve the node utilization but this is
      a project for another day...
    */
    void removeSubset(Subset const &subset);

    void print(bool decend = true) const;

    void getAllEntries(Subset &items,
                       Subset &duplicates);

    // alexr, 10-1-2013 moving debug utility here since we want to
    //                  use it all over the place
    void checkVsSubset(Subset const &ids);

    //uncomment this to turn on these functions for debugging
    //#define InsertionDebugOn
#ifdef InsertionDebugOn
    //chrisgdebug
    bool checkBoundingBoxes();
    bool checkLeafLevels();
    double calculateMargin();
    double calculateOverlap();
    void printLevelStats();
#endif

    class Data;

  protected:
    void checkForPrevQuery() const;

  private:

    Data *_data;

    // hide copy/assignment
    Rtree(Rtree const &rhs);
    Rtree &operator=(Rtree const &rhs);
  };
}

#endif
