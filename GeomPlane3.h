#ifndef _GEOM_PLANE3_H
#define _GEOM_PLANE3_H


#include "GeomVector.h"

namespace Geom
{
  class Line3;
  class Segment3;

  /*! The Plane3 is a geometry class storing a plane in 3-D. */
  class Plane3
  {
  public:
    Plane3()
      : _haseqn(false)
    {}

    /*! Construct a plane passing through \a point with normal \a normal.
      The normal vector MUST be unit length! */
    Plane3(Vector3d const &point,
           Vector3d const &normal)
      : _point(point)
      , _normal(normal)
      , _haseqn(false)
    {}

    /*! Construct a plane passing through three points.  The normal is oriented
     * in the direction p1-p2-p3 (right hand rule) */
    Plane3(Vector3d const &p1,
           Vector3d const &p2,
	   Vector3d const &p3)
      : _point(p1)
      , _normal()
      , _haseqn(false)
    {
      Vector3d vec1(p2-p1);
      Vector3d vec2(p3-p1);
      _normal = vec1.cross(vec2);
      _normal.normalize();
    }

    /*! Get the point on the plane. */
    Vector3d const & getPoint() const {return _point;}

    /*! Set the point on the plane. */
    void setPoint(Vector3d const &point) {
      _point = point;
      _haseqn = false;
    }

    /*! Get the normal of the plane. */
    Vector3d const &getNormal() const {return _normal;}
    /*! Set the normal on the plane.  The normal vector
      MUST be unit! */
    void setNormal(Vector3d const &normal) {
      _normal = normal;
      _haseqn = false;
      // fixme -- add development assertion that normal is unit-length
    }

    /*! Get the equation of the plane. */
    double const (&getEquation() const)[4] {
      if (!_haseqn)
        this->computeEquation();
      return _eqn;
    }
        
    double getConstant() const {
      if (!_haseqn)
        this->computeEquation();
      return _eqn[3];
    }

    /*! Find intersection with another plane in 3-D.  
      The result is a line (in 3-D).
    
      Returns false if the 2 planes are parallel, true otherwise.
    */
    bool findIntersection(Plane3 const &plane2,
			  Line3 &intLine) const;

    /*! Compute the signed distance from a plane to a point.  If the
      argument /a closestpoint is non-null, the closest point on the plane
      will be placed in the input vector. */
    double getSignedDistance(Vector3d const &point,
			     Vector3d *closestpoint = 0) const;

    /*! Compute the absolute distance from a plane to a point.  If the
      argument /a closestpoint is non-null, the closest point on the plane
      will be placed in the input vector. */
    inline double getDistance(Vector3d const &point,
			      Vector3d *closestpoint = 0) const
    {
      return abs(getSignedDistance(point, closestpoint));
    }

    /*!
      find intersection of a linesegment with a plane.
      if the segment lies completely in the plane returns true but segloc is -1
      otherwise segloc is always >=0  and <=1.0 when it returns true and is 
      the parametric distance the intersection occurs from coordinate1. 
      parametric distance means the intersection point can be computed
      by coord1 + segloc*(coord2-coord1) or seg.getCoordinate1() + 
      segloc*seg.getDirection() - this way the square root is avoided
      segloc is unchanged if it no intersection occurs.
    */
    //I updated this so segloc is the parametric distance rather than actual
    //distance as I explain above - this is consistent with all the other
    //functions and allows us to avoid taking a square root - chrisg 2/16/2005
    //NOTE:: if you want an exact precision version of plane/seg intersection
    // you can use the Triangle3::findExactIntersection with the infinitePlane
    //variable set to true (see Triangle3.h) This will then test the 
    //intersection of the infinite plane defined by the 3 points in the 
    //triangle with the segment rather than just the triangle
    //I have assigned myself a task to add an exact--or at least more exact
    //than this-- precision version of this. chrisg 8/9/2011
    bool findIntersection(Segment3 const &seg,
			  double* segloc) const;

    /*! Intersect plane with line in 3-D.  The result is the location along the
     *  line

     Returns false if the line is parallel to the plane, true otherwise.
    */
    bool findIntersection(Line3 const &line,
			  double* linloc) const;
  private:
    Vector3d _point;
    Vector3d _normal;
    mutable bool _haseqn;
    mutable double _eqn[4];

    void computeEquation() const;
  };


}

#endif
