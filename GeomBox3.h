#ifndef _Geom_Box3
#define _Geom_Box3


#include <cmath>
#include <vector>

#include "GeomBoundedShape3.h"

#include "GeomVector.h"

namespace Geom
{
  class BoundingBox;

  class ConeFrustum3;
  class Point3;
  class Ray3;
  class Segment3;
  class Sphere3;

/**
   * An axis-aligned box class in 3d.
   */
  class Box3 : public BoundedShape3
  {
  public:
    /**
     * Construct a box given two opposite corner points.  These points must
     * lie across a diagonal to properly specify the box but need not be the
     * actual min/max points.
     */
    explicit Box3(Vector3d const &corner1,
		  Vector3d const &corner2);

    /**
     * Construct a box given a set of points.
     */
    explicit Box3(std::vector < Vector3d > const &queryFace,
                  double adjustment = 0.0);

    /** Empty constructor.  Corner vectors will not be initialized. */
    Box3() {}

    Box3(BoundingBox const &bbox);
    
    /** Copy constructor. */
    Box3(Box3 const &rhs);


    /** Assignment operator. */
    Box3 &operator=(Box3 const &rhs)
    {
      _min = rhs._min;
      _max = rhs._max;
      return *this;
    }

    Shape3 *clone() const
    {
      return new Box3(*this);
    }

    /** Return the lower corner. */
    Vector3d const & min() const {return _min;}
    /** Return the lower corner. */
    Vector3d       & min()       {return _min;}

    /** Return the upper corner. */
    Vector3d const & max() const {return _max;}
    /** Return the upper corner. */
    Vector3d       & max()       {return _max;}

    /**
     * The position of the center of the box.
     */
    Vector3d getCenter() const {
      return Vector3d(0.5*(_min(0)+_max(0)),
                      0.5*(_min(1)+_max(1)),
                      0.5*(_min(2)+_max(2)));
    }
    
    /**
     * A vector from the minimum corner to the maximum corner.
     */
    Vector3d getDiagonal() const {
      return Vector3d(_max(0) - _min(0),
                      _max(1) - _min(1),
                      _max(2) - _min(2));
    }

    /**
     * The Extent (1/2 width) of the box from minimum to maximum.
     */
    double getExtent(int direction) const {
    return 0.5*(_max(direction) - _min(direction));
    }
    /**
     * Inflate the box in all directions, + and -, by the given amount so the 
     * resulting change in width is 2*inflateAmt
     */
    void inflate(double const inflateAmt);
    /**
     * Implement BoundedShape3::getBoundingBox() method.
     */
    BoundingBox getBoundingBox() const;

    /**
     * Implement BoundedShape3::intersectsBox() method.
     */
    bool intersectsBox(BoundingBox const &bbox) const;

    /**
     * Implement Shape::containsPoint() method.
     */
    bool containsPoint(Vector3d const &point) const;

    /**
     * Returns true if box1 and box2 overlap (with open ranges), false 
     * otherwise.
     */
    bool testIntersection(Box3 const &box2) const {
      if (this->min()(0) > box2.max()(0) ||
	  this->max()(0) < box2.min()(0))
	return false;
      if (this->min()(1) > box2.max()(1) ||
	  this->max()(1) < box2.min()(1))
	return false;
      if (this->min()(2) > box2.max()(2) ||
	  this->max()(2) < box2.min()(2))
	return false;
      return true;
    }

    /**
     * Returns true if the point is contained within the box (inclusive),
     * false otherwise.
     */
    bool testIntersection(Vector3d const &coord) const {
      return (coord(0) >= this->min()(0) &&
	      coord(1) >= this->min()(1) &&
	      coord(2) >= this->min()(2) &&
	      coord(0) <= this->max()(0) &&
	      coord(1) <= this->max()(1) &&
	      coord(2) <= this->max()(2));
    }
    
    /**
     * Return the distance squared between the box and the point.  If the
     * point is in the box, zero will be returned.
     */
    double getDistanceSquared(Vector3d const &coord) const;

    /**
     * Return the distance between the box and the point.  If the
     * point is in the box, zero will be returned.
     */
    inline double getDistance(Vector3d const &coord) const {
      return sqrt(getDistanceSquared(coord));
    }
    
    /**
     * Find the intersection between a Box3 and Segment3 object.  
     * If the segment pierces the box or is contained by it, true will be returned. 
     * If the segment does not intersect the box false will
     * be returned   
     *
     * This function could easily be modified to return the points of intersection
     * with the box as a distance from the first coordinate.
     *
     * NOTE:
     * This method has some tolerance issues if the segment runs through
     * a corner of the box.  As this method is used by octree for refinement
     * a segment might not intersect any box if it passes along a corner edge
     * along one axis.  so I have added an inflation to create a slightly 
     * inflated box before intesection to ensure that the segment intersects.  
     * inflateFactor is multiplied by the average dimension of the box
     * and the box is inflated by that amount -- by default 0.0.  When used
     * in the octree we use 0.001 for inflation.
     */
    bool testIntersection(Segment3 const &segment,
                          double inflateFactor=0.0) const;
  
    /**
     * Return the distance squared between the box and the segment. If the
     * segment is in the box, zero will be returned.
     */
    double getDistanceSquared(Segment3 const &segment3) const;

    /**
     * Return the distance between the box and the segment. If the segment is 
     * in the box, zero will be returned.
     */
    inline double getDistance(Segment3 const &segment3) const {
      return sqrt(getDistanceSquared(segment3));
    }  

    /**
     * Find the intersection between a Box3 and Ray3 object.  If the ray
     * pierces the box, true will be returned and the distance between the
     * ray origin and the closest intersection point will be returned in the 
     * distance argument.  If the ray does not pierce the box, false will be 
     * returned and the distance argument will be indeterminate.
     * NOTE:
     * This method has some tolerance issues if the ray runs through
     * a corner of the box.  along one axis.  
     * so I have added an inflation to create a slightly 
     * inflated box before intesection to ensure that the ray intersects.  
     * inflateFactor is multiplied by the average dimension of the box
     * and the box is inflated by that amount -- by default 0.0 no inflation
     */
    bool testIntersection(Ray3 const &ray,
			  double *distance,
                          double inflateFactor=0.0) const;
  
    /**
     * Returns true if the box intersects the sphere, false otherwise.
     */
    bool testIntersection(Sphere3 const &sphere) const;

    /**
     * Returns true if the box intersects the frustum, false otherwise.
     * WARNING!!!
     * This test is not completely precise.  It may get some bboxes that
     * are right around one end of the frustum but just outside the frustum
     * end cap.  This issue has been in here from the beginning so I am
     * not going to fix it the day before the cut 4.02 as I don't have a fix
     * worked out and it needs sufficient testing.  For the purposes of
     * refinement, it is fine but not precise and should be fixed.  
     *  I have assigned myself a bug report CCMP-6333
     * Chris Geisert 01/05/2008
     */
    bool testIntersection(ConeFrustum3 const &frustum) const;

    bool isDegenerate() const;

  private:    
    Vector3d _min;
    Vector3d _max;
  };

}

#endif
