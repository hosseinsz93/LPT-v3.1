#ifndef _GEOM_TRIANGLE3_H
#define _GEOM_TRIANGLE3_H

#include "GeomBoundedShape3.h"

#include "GeomPoint3.h"
#include "GeomSegment3.h"
#include "GeomLine3.h"
#include "GeomVector.h"

namespace Geom
{
  class Ray3;

  /*! The Triangle3 is a geometry class storing a triangle in 3-D. Things
    like the triangle area, normal, and unit normal are cached behind the 
    scenes so subsequent retrieval is a no-op. */
  class Triangle3 : public BoundedShape3
  {
  public:
    Triangle3()
      : BoundedShape3()
      ,	_hasedges(false)
      ,	_hasnormal(false)
      ,	_hasunitnormal(false)
      ,	_hasarea(false)
    {}

    Triangle3(Vector3d const (&v)[3])
      : BoundedShape3()
      , _hasedges(false)
      , _hasnormal(false)
      , _hasunitnormal(false)
      , _hasarea(false)
    {
      _vertices[0].position() = v[0];
      _vertices[1].position() = v[1];
      _vertices[2].position() = v[2];
    }

    Triangle3(Vector3d const &v1, 
              Vector3d const &v2,
              Vector3d const &v3)
      : BoundedShape3()
      , _hasedges(false)
      , _hasnormal(false)
      , _hasunitnormal(false)
      , _hasarea(false)
    {
      _vertices[0].position() = v1;
      _vertices[1].position() = v2;
      _vertices[2].position() = v3;
    }

    Triangle3(Triangle3 const &rhs)
      : BoundedShape3(rhs)
    {
      *this = rhs;
    }

    Triangle3 & operator=(Triangle3 const &rhs)
    {
      _vertices[0] = rhs._vertices[0];
      _vertices[1] = rhs._vertices[1];
      _vertices[2] = rhs._vertices[2];
      if ((_hasedges = rhs._hasedges)) { _edges[0] = rhs._edges[0]; _edges[1] = rhs._edges[1]; _edges[2] = rhs._edges[2]; }
      if ((_hasnormal = rhs._hasnormal)) { _normal = rhs._normal; }
      if ((_hasunitnormal = rhs._hasunitnormal)) { _unitnormal = rhs._unitnormal; }
      if ((_hasarea = rhs._hasarea)) { _area = rhs._area; }
      return *this;
    }

    Point3 const (& getPoints() const)[3] { return _vertices; }

    Point3 const & getPoint(int i) const { return getPoints()[i]; }

    Vector3d const & getCoordinate(int i) const { return getPoint(i).position(); }

    /** Set the coordinates of the triangle. */
    void setCoordinates(Vector3d const &v1, 
                        Vector3d const &v2,
                        Vector3d const &v3) {
      _vertices[0].position() = v1;
      _vertices[1].position() = v2;
      _vertices[2].position() = v3;
      _hasedges = false;      
      _hasnormal = false;
      _hasunitnormal = false;
      _hasarea = false;
    }

    /*! Get the n'th edge (0 <=n < 3) of the triangle, returned as a 
      Segment3 object. */
    Segment3 const & getEdge(int n) const {
      if (!_hasedges)
        computeEdges();
      return _edges[n];
    }

    void setIds(Id const &id1,
		Id const &id2,
		Id const &id3) {
      _vertices[0].setId(id1);
      _vertices[1].setId(id2);
      _vertices[2].setId(id3);
      if (_hasedges)
	{
	  _edges[0].setIds(id1,id2);
	  _edges[1].setIds(id2,id3);
	  _edges[2].setIds(id3,id1);
	}
    }

    void getIds(Id &id1,
		Id &id2,
		Id &id3) const {
      id1 = _vertices[0].getId();
      id2 = _vertices[1].getId();
      id3 = _vertices[2].getId();
    }

    /*! Get the (non-unit) normal. */
    Vector3d const &getNormal() const {
      if (!_hasnormal)
        this->computeNormal();
      return _normal;
    }
    /*! Set the (non-unit) normal. */
    void setNormal(Vector3d const &normal) {
      _normal = normal;
      _hasnormal = true;
    }
                  
    /*! Get the unit normal. */
    Vector3d const &getUnitNormal() const {
      if (!_hasunitnormal)
        this->computeUnitNormal();
      return _unitnormal;
    }
    /*! Set the unit normal. */
    void setUnitNormal(Vector3d const &unitnormal) {
      // fixme -- add development assertion that vector is unit length
      _unitnormal = unitnormal;
      _hasunitnormal = true;
    }

    /*! Get the triangle area. */
    double getArea() const {
      if (!_hasarea)
        this->computeUnitNormal(); // need it for area
      return _area;
    }
    /*! Set the triangle area. */
    void setArea(double area) {
      _area = area;
      _hasarea = true;
    }


    // Get triangle centroid -- currently computed every time!
    Vector3d getCentroid() const;


    /// Implement BoundedShape3 interface
    BoundingBox getBoundingBox() const;

    bool intersectsBox(BoundingBox const &bbox) const;

    /// Implement Shape3 virtual interface
    Shape3 *clone() const {
      return new Triangle3(*this);
    }

    bool containsPoint(Vector3d const &point) const;

    double getQuality() const;

    /*! Test intersection between a triangle and a segment in 3-D.  The 
      tolerance parameter "inflates" the triangle in all directions. Returns
      true if the triangle intersects the segment and false otherwise. */
    bool testIntersection(Segment3 const &segment,
			  double tolerance) const;

    /*! Test/find intersection using exact intersection algorithm.
      \todo these routines currently don't deal with planar intersections!
    */
    bool testExactIntersection(Segment3 const &segment,
			       bool infinitePlane=false) const;
    bool findExactIntersection(Segment3 const &segment,
			       Vector3d &IntPt,
			       bool infinitePlane=false) const;

    /*! Test/find intersection using exact intersection algorithm with simulation
      of simplicity.
    */
    bool testExactIntersection_sos(Segment3 const &segment,
				   bool infinitePlane=false) const;
    bool findExactIntersection_sos(Segment3 const &segment,
				   Vector3d &IntPt,
				   bool infinitePlane=false,
				   bool * segmentOrientation=0) const;

    /*! Test/find intersection using exact intersection algorithm.
      if the infinitePlane flag is set to true, the lhs triangle (this triangle) will be considered
      an infinite plane and the check will return whether or not triangle2 intersects that infinite plane
    */
    bool testExactIntersection(Triangle3 const &triangle2, bool infinitePlane = false) const;

    bool findExactIntersection(Triangle3 const &triangle2,
			       Segment3        &IntSeg) const;

    /*! Test/find intersection using exact intersection algorithm with simulation
      of simplicity.
    */
    bool testExactIntersection_sos(Triangle3 const &triangle2) const;

    bool findExactIntersection_sos(Triangle3 const &triangle2,
				   Segment3        &IntSeg) const;

    /*! Compute the squared distance from a triangle to a point. If the
      argument /a closestpoint is non-null, the closest point on the triangle
      will be placed in the input vector. */
    double getDistanceSquared(Vector3d const &point,
			      Vector3d *closestpoint = 0) const
      {
        return getDistanceSquared2(point,closestpoint);
      }
    
    double getDistanceSquared2(Vector3d const &point,
			       Vector3d *closestpoint = 0) const;
    
    /*! Compute the distance from a triangle to a point.  If the
      argument /a closestpoint is non-null, the closest point on the triangle
      will be placed in the input vector. */
    double getDistance(Vector3d const &point,
		       Vector3d *closestpoint = 0) const
    {
      return sqrt(getDistanceSquared(point, closestpoint));
    }

    /*! Compute the distance from a triangle to a line
     */
    double getDistanceSquared(Line3 const &line,
			      double *lineParam = 0,
			      Vector3d *closestTriPoint = 0) const;
  
    double getDistance(Line3 const &line,
		       double *lineParam = 0,
		       Vector3d *closestTriPoint = 0) const
    {
      return sqrt(getDistanceSquared(line,lineParam,closestTriPoint));
    }

    /*! Compute the distance from a triangle to a line setment. */
    double getDistanceSquared(Segment3 const &seg,
			      double *segParam = 0,
			      Vector3d *closestTriPoint = 0) const;

    double getDistance(Segment3 const &seg,
		       double *segParam = 0,
		       Vector3d *closestTriPoint = 0) const
    {
      return sqrt(getDistanceSquared(seg,segParam,closestTriPoint));
    }

    /*! Compute the distance between two triangles */
    double getDistanceSquared(Triangle3 const &tri,
			      Vector3d *closestTriPoint=0) const;
    double getDistance(Triangle3 const &tri) const
    {return sqrt(getDistanceSquared(tri));}


    /*! Return true if triangles overlap in plane of this triangle */
    bool overlaps(Triangle3 const &tri) const;


    /*! Return the height of the triangle (the distance from the longest edge
        to the opposite vertex)
    */
    double getHeightSquared(double *baseLenSquared=0) const;
    double getHeight(double *baseLen=0) const
    {
      double baseLenSq;
      double heightSq = this->getHeightSquared(&baseLenSq);
      if (baseLen)
	*baseLen = sqrt(baseLenSq);
      return sqrt(heightSq);
    }

    /*! return the radius of the largest inscribed circle
     */
    double getInradius() const;

    /*! Find intersection between a triangle and a ray in 3-D.  There is
      currently no tolerance parameter in this function -- do we need one?
    
      Returns true if the triangle intersects the ray and false otherwise.
    
      The distance parameter is an optional argument -- if non-zero, the 
      distance along the ray to the intersection point will be returned if
      the ray intersects the triangle.  No value will be stored if the ray
      does not intersect the triangle.
    */
    bool testIntersection(Ray3 const &ray,
			  double *distance = 0) const;
    bool testIntersection(Line3 const &line,
			  double *distance = 0) const; // signed distance!

    /*! Find intersection using exact intersection algorithm.  Optional \a Dot
        is the exact value of the dot product of the ray direction and the
        triangle normal.
    */
    bool testExactIntersection(Ray3 const &ray,
                               double *distance = 0,
                               double *Dot = 0) const;
    bool testExactIntersection(Line3 const &line,
                               double *distance = 0, // signed distance!
                               double *Dot = 0) const;

    /*! Solves for the barycentric coordinates of a vertex in the plane of the triangle.
     */
    void CartesianToBarycentricCoordinates(Vector3d const &p, double result[3]) const
    {
      return CartesianToBarycentricCoordinates(p(0),p(1),p(2),result);
    }
    
    void CartesianToBarycentricCoordinates(double pX, double pY, double pZ, double result[3]) const;

    void BarycentricToCartesianCoordinates(Vector3d const &b, double result[3]) const
    {
      return BarycentricToCartesianCoordinates(b(0),b(1),b(2),result);
    }
    void BarycentricToCartesianCoordinates(double bX, double bY, double bZ, double result[3]) const;

  private:
    Point3 _vertices[3];

    mutable bool _hasedges;
    mutable Segment3 _edges[3];

    mutable bool _hasnormal;
    mutable Vector3d _normal;
    mutable bool _hasunitnormal;
    mutable Vector3d _unitnormal;
    mutable bool _hasarea;
    mutable double _area;

    inline void computeEdges() const
    {
      _edges[0].setCoordinates(_vertices[0].position(), _vertices[1].position());
      _edges[1].setCoordinates(_vertices[1].position(), _vertices[2].position());
      _edges[2].setCoordinates(_vertices[2].position(), _vertices[0].position());
      _edges[0].setIds(_vertices[0].getId(), _vertices[1].getId());
      _edges[1].setIds(_vertices[1].getId(), _vertices[2].getId());
      _edges[2].setIds(_vertices[2].getId(), _vertices[0].getId());
      _hasedges = true;
    }

    inline void computeNormal() const
    {
      Vector3d const vec1(_vertices[2].position() - _vertices[1].position());
      Vector3d const vec2(_vertices[0].position() - _vertices[1].position());
      _normal = vec1.cross(vec2)*0.5;
      _hasnormal = true;
    }

    void computeUnitNormal() const;
  };

}

#endif
