#include "GeomPlane3.h"
#include "GeomPoint3.h"
#include "GeomRay3.h"
#include "GeomLine3.h"
#include "GeomTriangle3.h"

#include "GeomBoundingBox.h"
#include "GeomUtility.h"

#include "GeomFaceCalculations.h"

#include "Geompredicates.h"
#include "Geomorient3d_sos.h"

#include <limits>

namespace Geom
{
  static Line3 _lineBuffer;

  void Triangle3::computeUnitNormal() const
  {
    _unitnormal = getNormal();
    _area = _unitnormal.normalize();
    _hasunitnormal = true;
    _hasarea = true;
  }

  Vector3d Triangle3::getCentroid() const
  {
    return triangleCentroid(this->getCoordinate(0),
			    this->getCoordinate(1),
			    this->getCoordinate(2));
  }

  BoundingBox Triangle3::getBoundingBox() const
  {
    BoundingBox bbox(_vertices[0].position());
    bbox.expand(_vertices[1].position());
    bbox.expand(_vertices[2].position());
    return bbox;
  }

  bool Triangle3::intersectsBox(BoundingBox const &bbox) const
  {
    // consider migrating this to a call in Box3Triangle3
    return bbox.intersectsTriangle(_vertices[0].position(),
				   _vertices[1].position(),
				   _vertices[2].position());
  }

  bool Triangle3::containsPoint(Vector3d const &point) const
  {
    return (getDistanceSquared(point) <= ALMOST_ZERO);
  }


#if 1 // this is the new quality measure
  static inline double perimeter(Triangle3 const &Triangle)
  {
    return (Triangle.getEdge(0).getLength() +
            Triangle.getEdge(1).getLength() +
            Triangle.getEdge(2).getLength());
  }

  static inline double halfperimeter(Triangle3 const &Triangle)
  {
    return 0.5 * perimeter(Triangle);
  }

  static inline double inradius(Triangle3 const &Triangle)
  {
    // Compute the triangle inradius
    double halfperim = halfperimeter(Triangle);
    if (halfperim < ALMOST_ZERO)
      return 0.0;

    return Triangle.getArea()/halfperim;
  }

  double Triangle3::getInradius() const
  {
    return inradius(*this);
  }

  static inline double circumradius(Triangle3 const &Triangle)
  {
    // Compute the triangle circumradius
    double area = Triangle.getArea();
    if (area < ALMOST_ZERO)
      return std::numeric_limits<double>::max();

    return (Triangle.getEdge(0).getLength() *
	    Triangle.getEdge(1).getLength() *
	    Triangle.getEdge(2).getLength() /
	    (4.0*area));
  }

  double Triangle3::getQuality() const
  {
    /* Compute the quality of this triangle.  The quality is given by:
     * Quality = factor*(r/R)
     * where:
     * factor = 2.0
     * r = triangle inradius = A/s
     * R = triangle circumradius = abc/(4A)
     * A = triangle area
     * s = triangle half perimeter
     */
    double den = circumradius(*this);

    if (den < ALMOST_ZERO)
      return 0.0;

    return 2.0*inradius(*this)/den;
  }

#else // this is the old quality measure
  double Triangle3::getQuality() const
  {
    /* Compute the quality of this triangle.  The quality is given by:
     * Quality = factor*A/(R*R)
     * where:
     * factor = 4*sqrt(3)/9 = 0.769800358
     * A = triangle area
     * R = triangle circumradius
     * R = abc/(4A)
     * a,b,c = triangle edge lengths
     * Simplifying, we get:
     * Quality = factor2 * A*A*A/(a*a*b*b*c*c)
     * where:
     * factor2 = 16*factor
     * So, we basically need all the lengths squared and the area...
     */

    double const factor2 = 12.31680574;

    double areamag = getArea();

    // Return 0 quality if 0 area triangle
    if (areamag < ALMOST_ZERO)return 0.0;

    Vector3d edge1 = getEdge1().getDirection();
    Vector3d edge2 = getEdge2().getDirection();
    Vector3d edge3 = getEdge3().getDirection();

    double elen1sq = edge1.dot(edge1);
    double elen2sq = edge2.dot(edge2);
    double elen3sq = edge3.dot(edge3);

    double num = factor2*areamag*areamag*areamag;
    double den = elen1sq*elen2sq*elen3sq;

    return num/den;
  }
#endif


  bool Triangle3::testExactIntersection(Triangle3 const &triangle2, bool infinitePlane) const
  {
    // Test segments of triangle2 with triangle1
    if (testExactIntersection(triangle2.getEdge(0), infinitePlane)) return true;
    if (testExactIntersection(triangle2.getEdge(1), infinitePlane)) return true;
    if (testExactIntersection(triangle2.getEdge(2), infinitePlane)) return true;

    //only continue the test if we did not check for intersection on an infinite plane
    if (infinitePlane)
      return false;

    // Test segments of triangle1 with triangle2
    if (triangle2.testExactIntersection(getEdge(0))) return true;
    if (triangle2.testExactIntersection(getEdge(1))) return true;
    if (triangle2.testExactIntersection(getEdge(2))) return true;

    return false;
  }




  bool Triangle3::findExactIntersection(Triangle3 const &triangle2,
					Segment3 &IntSeg) const
  {
    unsigned n=0;
    Vector3d p[6], pAux;

    // Segments of this with triangle 2
    if (triangle2.findExactIntersection(getEdge(0), pAux))
      {p[n] = pAux; ++n;}					 
    if (triangle2.findExactIntersection(getEdge(1), pAux))
      {p[n] = pAux; ++n;}						 
    if (triangle2.findExactIntersection(getEdge(2), pAux))
      {p[n] = pAux; ++n;}

    // Segments of triangle 2 with this
    if (findExactIntersection(triangle2.getEdge(0), pAux))
      {p[n] = pAux; ++n;}						 
    if (findExactIntersection(triangle2.getEdge(1), pAux))
      {p[n] = pAux; ++n;}						 
    if (findExactIntersection(triangle2.getEdge(2), pAux))
      {p[n] = pAux; ++n;}

    // TODO: We should consider what to do about n > 2
    if (n < 2)
      return false;
    else
      IntSeg.setCoordinates(p[0],p[1]);
    return true;
  }




  bool Triangle3::testExactIntersection_sos(Triangle3 const &triangle2) const
  {
    // Test segments of this with triangle2
    if (triangle2.testExactIntersection_sos(getEdge(0))) return true;
    if (triangle2.testExactIntersection_sos(getEdge(1))) return true;
    if (triangle2.testExactIntersection_sos(getEdge(2))) return true;

    // Test segments of triangle2 with this
    if (testExactIntersection_sos(triangle2.getEdge(0))) return true;
    if (testExactIntersection_sos(triangle2.getEdge(1))) return true;
    if (testExactIntersection_sos(triangle2.getEdge(2))) return true;

    return false;
  }




  bool Triangle3::findExactIntersection_sos(Triangle3 const &triangle2,
					    Segment3 &IntSeg) const
  {
    unsigned n=0;
    Vector3d p[6], pAux;

    // Segments of triangle 1 with triangle 2
    if (triangle2.findExactIntersection_sos(getEdge(0), pAux))
      {p[n] = pAux; ++n;}						 
    if (triangle2.findExactIntersection_sos(getEdge(1), pAux))
      {p[n] = pAux; ++n;}						 
    if (triangle2.findExactIntersection_sos(getEdge(2), pAux))
      {p[n] = pAux; ++n;}

    // Segments of triangle 2 with triangle 1
    if (findExactIntersection_sos(triangle2.getEdge(0), pAux))
      {p[n] = pAux; ++n;}						 
    if (findExactIntersection_sos(triangle2.getEdge(1), pAux))
      {p[n] = pAux; ++n;}						 
    if (findExactIntersection_sos(triangle2.getEdge(2), pAux))
      {p[n] = pAux; ++n;}

    // TODO: We should consider whether it's possible to have n > 2 and how to handle that case
    if (n != 2)return false;

    IntSeg.setCoordinates(p[0],p[1]);
    return true;
  }
#if 0
//I am going to set this in the header file to call getDistanceSquared2
//Aly wrote getDistanceSquared2 because he saw problems in the
//trianglePointDistanceSquared -- now I am seeing problems with long
//skinny triangles so with Aly's agreement, I am changing all usage
//to getDistanceSquared2.
//If it no problems occur we can rename getDistanceSquared2 to 
//getDistanceSquared -- (geisert 09/16/2009)
  double Triangle3::getDistanceSquared(Vector3d const &querypt,
				       Vector3d *closestpoint) const
  {

    return trianglePointDistanceSquared(getCoordinate(0),
					getCoordinate(1),
					getCoordinate(2),
					querypt,
					closestpoint);

  }
#endif

  /* I (Aly) had precision problems with the above routine when a point was
   * on a triangle (returned non-zero distance).  I have rewritten it
   * my way and it seems to work.  I will just call it something different and
   * if it seems to be okay for a while, we can just remove the above routine.
   * 3/27/2005
   */
  double Triangle3::getDistanceSquared2(Vector3d const &point,
					Vector3d *closestpoint) const
  {
    Vector3d normal=getNormal();
    Vector3d cross =
      normal.cross(getEdge(0).getDirection());
    Vector3d vec(point - getCoordinate(0));
    if (cross.dot(vec) > 0.0)
      {
	cross = normal.cross(getEdge(1).getDirection());
	vec   = point - getCoordinate(1);
	if (cross.dot(vec) > 0.0)
	  {
	    cross=normal.cross(getEdge(2).getDirection());
	    vec   = point - getCoordinate(2);
	    if (cross.dot(vec) > 0.0)
	      {
		// Closest point is inside triangle
		double dot = normal.dot(vec);
		double magsq = normal.dot(normal);
		if (closestpoint)
		  {
		    double ratio = dot/magsq;
		    *closestpoint = point - normal*ratio;
		  }
		return dot*dot/magsq;
	      }
	  }
      }

    // Closest point is outside triangle.  Must be closest to one of the
    // segments of the triangle.  I can probably be smart about which ones to
    // check by retaining the above information but oh well...
    double dminsq = getEdge(0).getDistanceSquared(point,closestpoint);
    Vector3d loc;
    double dsq = getEdge(1).getDistanceSquared(point,&loc);
    if (dsq < dminsq)
      {
	dminsq = dsq;
	if (closestpoint)
	  *closestpoint = loc;
      }
    dsq = getEdge(2).getDistanceSquared(point,&loc);
    if (dsq < dminsq)
      {
	dminsq = dsq;
	if (closestpoint)
	  *closestpoint = loc;
      }
    return dminsq;
  }

  
  
  void Triangle3::
  CartesianToBarycentricCoordinates(double pX, double pY, double pZ, double result[]) const
  {
    double a,b,c,d,e,f,g,h,i;

    Vector3d v1, v2, v3;
    v1 = getPoint(0).position();
    v2 = getPoint(1).position();
    v3 = getPoint(2).position();

    a = v1(0) - v3(0);
    b = v2(0) - v3(0);
    c = v3(0) - pX;
    
    d = v1(1) - v3(1);
    e = v2(1) - v3(1);
    f = v3(1) - pY;
    
    g = v1(2) - v3(2);
    h = v2(2) - v3(2);
    i = v3(2) - pZ;
    
    /*
      Explanation for the following conditional statement:

      FIRST LINE: This is a necessary fix for when the triangle is perpendicular
      to the x-axis. We solve the problem by rotating the shape and point. Failure
      to use this may/will result in floating point exceptions being thrown

      SECOND AND THIRD LINE: These lines are an effort to reduce rounding error
      associated with having a small denominator when calculating result[0] and result[1]
      below. The denominator is calculated with the problem "rotated" to align with two 
      different axes, and if the "rotated" situation results in a larger denominator, 
      the problem is rotated. NB: This may introduce a significant trade off of speed for 
      accuracy.
     */

    if ((fabs(a) < ALMOST_ZERO && fabs(b) < ALMOST_ZERO) ||
	((fabs(d*(b+h) - e*(a+g)) > fabs(a*(e+h) - b*(d+g))) &&
	 (fabs(e*(a+g) - d*(b+h)) > fabs(b*(d+g) - a*(e+h)))))
      {
	double tmp = a; 
	a = d;
	d = tmp;
	
	tmp = b;
	b = e;
	e = tmp;
	
	tmp = c;
	c = f;
	f = tmp;	

      }
    
    result[0] = ((b*(f+i) - c*(e+h))/(a*(e+h) - b*(d+g)));
    result[1] = ((a*(f+i) - c*(d+g))/(b*(d+g) - a*(e+h)));
    result[2] = 1 - result[0] - result[1];

  }
  
  void Triangle3::
  BarycentricToCartesianCoordinates(double bX, double bY, double bZ, double result[]) const
  {
    Vector3d v1,v2,v3;
    v1 = getPoint(0).position();
    v2 = getPoint(1).position();
    v3 = getPoint(2).position();

    for (int i=0;i<3;i++)
      result[i] = bX*v1(i) + bY*v2(i) + bZ*v3(i);
  }


  double Triangle3::
  getDistanceSquared(Line3 const &line,
		     double* lineParam,
		     Vector3d *closestTriPoint) const
    
  {
    double lineDist;
    if (this->testIntersection(line,&lineDist)) // line intersects triangle
      {
        double mag = line.getDirection().mag();
        double tLineParam = lineDist/mag;
	Vector3d intpt(line.getOrigin() + 
			       Vector3d(line.getDirection()*tLineParam));

        if (lineParam) *lineParam = tLineParam;
        if (closestTriPoint) *closestTriPoint = intpt;

        return 0.0;//intersects so distance must be zero
      }

    // Now, we are either parallel, or the line intersects the plane of the
    // triangle, but not the triangle.  In either case, the closest distance
    // is the smallest of the distances to the edges of the triangle
	
    //initialize to prevent compiler warnings
    double closestLineParam=0.0,closestSegParam=0.0;
    int closestSegIndex=0;

    double minTriangleLineDistance = std::numeric_limits<double>::max(); 
	
    for (int i=0;i<3;i++) {
      double tSegParam,tLineParam;
      Segment3 const &triSeg = this->getEdge(i);
      double thisDistanceSqr = line.getDistanceSquared(triSeg,
                                                       &tLineParam,&tSegParam);
      if (thisDistanceSqr < minTriangleLineDistance) 
      {
        minTriangleLineDistance = thisDistanceSqr;
        if (lineParam)
          closestLineParam = tLineParam;
	  
        if (closestTriPoint)
        {
          closestSegIndex = i;
          closestSegParam = tSegParam;
        }
      }
    }
	
    if (lineParam) *lineParam = closestLineParam;
    if (closestTriPoint)
    {
      Segment3 const &closeSeg=this->getEdge(closestSegIndex);
      *closestTriPoint = 
        closeSeg.getCoordinate(0) + closeSeg.getDirection()*closestSegParam;
    }
    return minTriangleLineDistance;
  }

  double Triangle3::
  getDistanceSquared(Segment3 const &seg,
		     double *segParam,
		     Vector3d *closestTriPoint) const
  {
    Vector3d closeTriPoint, closeSegPoint;
    

    
    Vector3d dir = seg.getDirection();
    //if segment is zero length just use endpoint distance
    if (dir.mag2() < ALMOST_ZERO)
    {
      if (segParam)
        *segParam=0.0;
      return getDistanceSquared(seg.getPoint(0).position(),closestTriPoint);
    }

    double closeLineParam =0;
    double fSqrDist = getDistanceSquared(Line3(seg.getPoint(0).position(),dir),
					 &closeLineParam,
					 &closeTriPoint);
    
    if (closeLineParam < 0)
      { //triangle is closest to point(0) on segment
	fSqrDist = getDistanceSquared(seg.getPoint(0).position(),closestTriPoint);
	if (segParam) *segParam = 0;
      }
    else if ( closeLineParam > 1)
      { //triangle is closest to point(1) on segment
	fSqrDist = getDistanceSquared(seg.getPoint(1).position(),closestTriPoint);
	if (segParam) *segParam = 1;
      }
    else
      {// here, the closest point found on the line containing the segment is somewhere on the segment
	if (segParam) *segParam = closeLineParam;
	if (closestTriPoint) *closestTriPoint = closeTriPoint;
      }
    
    return fSqrDist;
  }




  double Triangle3::getDistanceSquared(Triangle3 const &tri,
				       Vector3d * closestTriPoint) const
  {
    // The distance between the two triangles is the shortest of the triangle
    // segment distances (done both ways)
    Vector3d * tmpClosestPointPtr = 0;
    Vector3d tmpClosestPoint;

    if (closestTriPoint)
      tmpClosestPointPtr = &tmpClosestPoint;

    double distSq = this->getDistanceSquared(tri.getEdge(0),NULL,tmpClosestPointPtr);
    if (closestTriPoint)
      *closestTriPoint = tmpClosestPoint;

    double temp = this->getDistanceSquared(tri.getEdge(1),NULL,tmpClosestPointPtr);
    if (temp < distSq)
      {
	distSq = temp;
	if (closestTriPoint)
	  *closestTriPoint = tmpClosestPoint;
      }
    temp = this->getDistanceSquared(tri.getEdge(2),NULL,tmpClosestPointPtr);
    if (temp < distSq)
      {
	distSq = temp;
	if (closestTriPoint)
	  *closestTriPoint = tmpClosestPoint;
      }

    temp = tri.getDistanceSquared(this->getEdge(0),NULL,tmpClosestPointPtr);;
    if (temp < distSq)
      {
	distSq = temp;
	if (closestTriPoint)
	  *closestTriPoint = tmpClosestPoint;
      }

    temp = tri.getDistanceSquared(this->getEdge(1),NULL,tmpClosestPointPtr);
    if (temp < distSq)
      {
	distSq = temp;
	if (closestTriPoint)
	  *closestTriPoint = tmpClosestPoint;
      }

    temp = tri.getDistanceSquared(this->getEdge(2),NULL,tmpClosestPointPtr);
    if (temp < distSq)
      {
	distSq = temp;
	if (closestTriPoint)
	  *closestTriPoint = tmpClosestPoint;
      }


    return distSq;
  }




  bool Triangle3::overlaps(Triangle3 const &tri) const
  {
    /* Does tri overlap this triangle in plane of this triangle?
     * For this to happen, there must be at least one vertex of tri in the
     * right half plane of each segment of this triangle
     */
    for (unsigned i=0; i<3; ++i) // loop over segments of this triangle
      {
	Segment3 const &seg = this->getEdge(i);
	Vector3d plane = this->getNormal().cross(seg.getDirection());
	bool found = false;
	for (unsigned j=0; j<3; ++j) // loop over verts of tri
	  {
	    Vector3d vec(tri.getCoordinate(j) - seg.getCoordinate(0));
	    if (plane.dot(vec) > 0) // note -- not enough to just touch!
	      {
		found = true;
		break;
	      }
	  }
	if (!found)
	  return false;
      }
    return true;
  }





  double Triangle3::getHeightSquared(double *baseLenSquared) const
  {
    unsigned longSeg = 0;
    double longSegLenSq = 0.0;
    for (unsigned i=0; i<3; ++i)
      {
	double lenSq = this->getEdge(i).getLengthSquared();
	if (lenSq > longSegLenSq)
	  {
	    longSegLenSq = lenSq;
	    longSeg = i;
	  }
      }

    if (baseLenSquared)
      *baseLenSquared = longSegLenSq;

    unsigned v = ((longSeg + 2) % 3);
    return this->getEdge(longSeg).getDistanceSquared(this->getCoordinate(v));
  }







  // local function needed for Segment->Triangle intersection
  static bool halfPlaneTest(Segment3 const &s,
                            Vector3d const &point,
                            Vector3d const &normal,
                            double tol,
                            double *kmin,
                            double *kmax)
  {
    Vector3d vec(s.getCoordinate(0) - point);
    double dot1 = normal.dot(vec);
    bool side1 = (dot1 >= -tol);

    vec = s.getCoordinate(1) - point;
    double dot2 = normal.dot(vec);
    bool side2 = (dot2 >= -tol);

    if (!side1 && !side2)
      return false;
    else if (side1 && side2)
      {
        *kmin = 0.0;
        *kmax = 1.0;
      }
    else if (!side1)
      {
        *kmax = 1.0;
        *kmin = (-tol - dot1) / (dot2 - dot1);
      }
    else 
      {
        *kmin = 0.0;
        *kmax = (-tol - dot1) / (dot2 - dot1);
      }
    return true;
  }


  // Another local function needed for Segment->Triangle intersection
  static bool boundsCheck(Segment3 const &s,
                          Vector3d const &p,
                          Vector3d const &vec,
                          double xmin, double xmax, double tol,
                          double *kmin, double *kmax)
  {
    /* Check what part (if any) of segment 's' lies inside the line defined by
     * p and vec with bounds xmin and xmax within tolerance tol.  Return the
     * parameter range that satisfies the bounds
     */

    double minval = xmin - tol;
    double maxval = xmax + tol;

    // Determine projection of segment on line
    Vector3d vec1(s.getCoordinate(0) - p);
    double     dot1 = vec.dot(vec1);

    Vector3d vec2(s.getCoordinate(1) - p);
    double     dot2 = vec.dot(vec2);

    // Determine sidedness with respect to minval
    bool s1minside = (dot1 >= minval);
    bool s2minside = (dot2 >= minval);
    if (!s1minside && !s2minside)return false; // both on -ve side

    // Determine sidedness with respect to maxval
    bool s1maxside = (dot1 <= maxval);
    bool s2maxside = (dot2 <= maxval);
    if (!s1maxside && !s2maxside)return false; // both on +ve side

    // At this point there is some overlap
    if (s1minside && s1maxside && s2minside && s2maxside) // both ends in bounds
      {
        *kmin = 0.0;
        *kmax = 1.0;
        return true;
      }

    // Trim the ends of the segment as needed (and store parameter values)
    if (!s1minside)
      *kmin = (minval - dot1) / (dot2 - dot1);
    else if (!s1maxside)
      *kmin = (maxval - dot1) / (dot2 - dot1);
    else
      *kmin = 0.0;

    if (!s2minside)
      *kmax = (minval - dot1) / (dot2 - dot1);
    else if (!s2maxside)
      *kmax = (maxval - dot1) / (dot2 - dot1);
    else
      *kmax = 1.0;

    return true;
  }



  bool Triangle3::testIntersection(Segment3 const &segment,
				   double tolerance) const
  {
    // use intersection routine from ammbatch aly_intersection.cpp
        
    /* Intersect triangle Triangle with edge Edge to within the specified
     * tolerance tol.  Return the intersection parameter values along the edge
     * (kmin and kmax).  Return value of function:
     *
     * DOES_NOT_INTERSECT  - no intersection between triangle and edge
     * RANGE_INTERSECTION  - edge intersects triangle from kmin to kmax
     *
     * NOTE:  This routine uses dimensional reduction to solve the problem
     *        exactly to within the user specified problem.  As a result, it
     *        never returns unique intersection, but always a range.
     *
     * Approach:
     * 1) 3-D problem:
     *    Does the edge intersect the "tolerance band" of the plane defined by
     *    the triangle?  If not, no intersection.  If so, then we retain the
     *    part of the edge that does and project it to the plane of the triangle.
     * 2) 2-D problem:
     *    Does any part of the edge lie on the right half of all three triangle
     *    edges (buffered with the tolerance band)?  If not, no intersection.
     *    If so, we retain the part of the edge that does and project it to each
     *    edge to ensure 1-D intersection is satisfied too.
     * 3) 1-D problem:
     *    Does the edge intersect the triangle edge buffered with the user
     *    specified tolerance?  If so, there is intersection.  Keep in mind here
     *    that if the triangle is obtuse, the intersection is checked with the
     *    extrema along that triangle edge.  This should be done by the calling
     *    program.
     */
        
    // Only continue if the triangle has non-zero area
    if (getArea() <= ALMOST_ZERO)
      return false;


    /* Get signed distances of edge points from triangle plane
     * If you are not within the tolerance band, you cannot have intersection...
     * e*_offplane =  0 -- on plane
     * e*_offplane = -1 -- on negative side of plane
     * e*_offplane = +1 -- on positive side of plane
     */
    Plane3 plane(getCoordinate(0), getUnitNormal());
    double d0 = plane.getSignedDistance(segment.getCoordinate(0));
    double d1 = plane.getSignedDistance(segment.getCoordinate(1));
    int e0_offplane = 0;
    if (d0 > tolerance)
      e0_offplane++;
    else if (d0 < -tolerance)
      e0_offplane--;

    int e1_offplane = 0;
    if (d1 > tolerance)
      e1_offplane++;
    else if (d1 < -tolerance)
      e1_offplane--;

    if (e0_offplane*e1_offplane > 0) return false;

    /* We now know that the edge segment intersects the tolerance band of the
     * plane that is defined by the triangle.  Clip the unwanted parts off and
     * project the remaining edge onto the plane of the triangle and perform 2D
     * intersection...
     *
     * NOTE:  No divide by 0 problem below because the only way d1 = d0 is if
     *        both points are inside the tolerance band, in which case the
     *        division is not performed...
     */

    double sd0 = d0;
    if (e0_offplane > 0)
      sd0 = tolerance;
    else if (e0_offplane < 0)
      sd0 = -tolerance;

    double k0;
    if (e0_offplane)
      k0 = (sd0 - d0)  / (d1 - d0);
    else
      k0 = 0.0;

    double sd1 = d1;
    if (e1_offplane > 0)
      sd1 = tolerance;
    else if (e1_offplane < 0)
      sd1 = -tolerance;
        
    double k1;
    if (e1_offplane)
      k1 = (sd1 - d0)  / (d1 - d0);
    else
      k1 = 1.0;

    Vector3d s0  = segment.getCoordinate(0);
    s0 += segment.getDirection()*k0;
    s0 -= this->getUnitNormal()*sd0;

    Vector3d s1 = segment.getCoordinate(0);
    s1 += segment.getDirection()*k1;
    s1 -= this->getUnitNormal()*sd1;

    Segment3 projected_edge(s0, s1);

    // projected_edge is now in plane of triangle, proceed to 2d check
    /*
     * Note that the edge can only intersect if it is on the right half of all
     * three edges of the triangle...
     */

    // Do half plane test for edge 0 of triangle
    Vector3d normal = this->getUnitNormal().cross(this->getEdge(0).getDirection());
    normal.normalize();
    double k0min, k0max;            
    if (!halfPlaneTest(projected_edge, this->getCoordinate(1), 
                       normal, tolerance, &k0min, &k0max))
      return false;

    normal = 
      this->getUnitNormal().cross(this->getEdge(1).getDirection());
    normal.normalize();
    double k1min, k1max;            
    if (!halfPlaneTest(projected_edge, this->getCoordinate(2), 
                       normal, tolerance, &k1min, &k1max))
      return false;               
        
    normal = 
      this->getUnitNormal().cross(this->getEdge(2).getDirection());
    normal.normalize();
    double k2min, k2max;            
    if (!halfPlaneTest(projected_edge, this->getCoordinate(0), 
                       normal, tolerance, &k2min, &k2max))
      return false;               
        
        
    // The portion of the edge that intersects the triangle is the 
    // intersection of all 3 ranges...

    double kmin, kmax;
    if (k0min > k1min)
      kmin = k0min;
    else
      kmin = k1min;
    if (kmin < k2min)
      kmin = k2min;
        
    if (k0max < k1max)
      kmax = k0max;
    else
      kmax = k1max;
    if (kmax > k2max)
      kmax = k2max;
        
    if (kmax < kmin)
      return false;


    /* We need to also remove the "tip-effects" (imagine a triangle with one
     * small angle -- adding tolerances to both edges attached to this angle
     * causes stuff to be found very far away from the vertex...), for which we
     * need to pass info to a 1-D intersection routine.
     */

    double xmin, xmax;
    double dot;
    Vector3d vec;

    // Bound in the direction of e1
    xmin = 0.0;
    xmax = this->getEdge(0).getLength();
    vec = this->getCoordinate(2) - this->getCoordinate(0);
    dot = this->getEdge(0).getUnitDirection().dot(vec);
    if (dot < xmin)
      xmin = dot;
    else if (dot > xmax)
      xmax = dot;
    if (!boundsCheck(projected_edge, this->getCoordinate(0),
                    this->getEdge(0).getUnitDirection(),
                    xmin, xmax, tolerance, &k0min, &k0max)) return false;

    if (k0min > kmin)
      kmin = k0min;
    if (k0max < kmax)
      kmax = k0max;

    // Bound in the direction of e2
    xmin = 0.0;
    xmax = this->getEdge(1).getLength();
    vec = this->getCoordinate(0) - this->getCoordinate(1);
    dot = this->getEdge(1).getUnitDirection().dot(vec);
    if (dot < xmin)
      xmin = dot;
    else if (dot > xmax)
      xmax = dot;
    if (!boundsCheck(projected_edge, this->getCoordinate(1),
                    this->getEdge(1).getUnitDirection(),
                    xmin, xmax, tolerance, &k1min, &k1max)) return false;

    if (k1min > kmin)
      kmin = k1min;
    if (k1max < kmax)
      kmax = k1max;

    // Bound in the direction of e3
    xmin = 0.0;
    xmax = this->getEdge(2).getLength();
    vec = this->getCoordinate(1) - this->getCoordinate(2);
    dot = this->getEdge(2).getUnitDirection().dot(vec);
    if (dot < xmin)
      xmin = dot;
    else if (dot > xmax)
      xmax = dot;
    if (!boundsCheck(projected_edge, this->getCoordinate(2),
                    this->getEdge(2).getUnitDirection(),
                    xmin, xmax, tolerance, &k2min, &k2max)) return false;

    if (k2min > kmin)
      kmin = k2min;
    if (k2max < kmax)
      kmax = k2max;
        
    if (kmax < kmin)
      {
        //      std::cout << "TIP EFFECTS!" << std::endl;
        return false;
      }

    return true;
  }




  // Exact intersection routines

  //
  // alexr, 1-24-2014, CCMP-59470
  //     Updating this routine to check each cartesian direction
  //     if the direction selected first fails since there are cases
  //     where the original routine was failing to find a point out
  //     of the plane, [e.g, (0,0,0),(0,0,.5),(1,1,0)].
  //     If there is any further trouble with this routine, we
  //     probably just ought to switch to using the normal to find
  //     the out of plane point, but we will keep it this way to
  //     avoid any general behavior changes or new robustness issues.
  //
  static bool findOutOfPlanePoint(Vector3d const &t0,
                                  Vector3d const &t1,
                                  Vector3d const &t2,
                                  Vector3d const &s0,
                                  Vector3d const &s1,
                                  Vector3d &P)
  {
    // The best out of plane point is a point in the direction that has the
    // minimum deviation (smallest bounding box cartesiand direction)
    // [alexr, 1-24-2014, this is not always true... see comment above this routine]
    BoundingBox bbox(t0);
    bbox.expand(t1);
    bbox.expand(t2);
    bbox.expand(s0);
    bbox.expand(s1);

    static Vector3d perturbX(1,0,0);
    static Vector3d perturbY(0,1,0);
    static Vector3d perturbZ(0,0,1);

    double len0 = bbox.sideLength(0);
    double len1 = bbox.sideLength(1);
    double len2 = bbox.sideLength(2);

    Vector3d *perturb0 = &perturbX;
    Vector3d *perturb1 = &perturbY;
    Vector3d *perturb2 = &perturbZ;

    // sort the three lengths
    if (len0 > len1)
      {
        std::swap(perturb0, perturb1);
        std::swap(len0, len1);
      }
    if (len1 > len2)
      {
        std::swap(perturb1, perturb2);
        std::swap(len1, len2);
        if (len0 > len1)
          {
            std::swap(perturb0, perturb1);
            std::swap(len0, len1);
          }
      }

    // Make sure you find a point P on the positive side of the triangle
    P = t0 + *perturb0;
    double len = orient3d((double*) &t0, (double*) &t1, (double*) &t2, (double*) &P);
    if (len > 0)
      return true;
    else if (len < 0)
      {
        P = t0 - *perturb0;
        return true;
      }
    else
      {
        // if first direction didn't get out of the plane, try another one
        P = t0 + *perturb1;
        len = orient3d((double*) &t0, (double*) &t1, (double*) &t2, (double*) &P);
        if (len > 0)
          return true;
        else if (len < 0)
          {
            P = t0 - *perturb1;
            return true;
          }
        else
          {
            // if we haven't succeeded yet, try the last direction
            P = t0 + *perturb2;
            len = orient3d((double*) &t0, (double*) &t1, (double*) &t2, (double*) &P);
            if (len > 0)
              return true;
            else if (len < 0)
              {
                P = t0 - *perturb2;
                return true;
              }
          }
      }
    // none of the directions worked so our triangle must be colinear...
    return false;

  }

  static bool findOutOfPlanePoints(Vector3d const &t0,
                                   Vector3d const &t1,
                                   Vector3d const &t2,
                                   Vector3d const &s0,
                                   Vector3d const &s1,
                                   Vector3d &P0,
                                   Vector3d &P1)
  {
    if (findOutOfPlanePoint(t0,t1,t2,s0,s1,P0))
      {
        P1 = t0 + (t0 - P0);
        return true;
      }
    return false;
  }

  static bool findOutOfPlanePointsUsingNormal(
      Vector3d const &t0,
      Vector3d const &t1,
      Vector3d const &t2,
      Vector3d &P0,
      Vector3d &P1)
  {
    Vector3d v1(t1-t0);
    Vector3d v2(t2-t0);
    Vector3d normal = v1.cross(v2);
    if (normal.mag2() < ALMOST_ZERO)
      return false;
    normal.normalize();

    Vector3d centroid =  triangleCentroid(t0,t1,t2);

    P0 = centroid + normal;
    P1 = centroid - normal;

    return true;
  }

  static bool exactIntersection(Vector3d const &t0,
				Vector3d const &t1,
				Vector3d const &t2,
				Vector3d const &s0,
				Vector3d const &s1,
				bool testCoplanar,
				Vector3d *IntPt=0,
				bool infinitePlane=false);

  static bool testEdge(Vector3d const &s0,
		       Vector3d const &s1,
		       Vector3d const &P,
		       Vector3d const &tA,
		       Vector3d const &tB,
		       Vector3d *IntPt,
		       std::vector< Vector3d > &Pts)
  {
    // Form a triangle with the segment and the out of plane point.  Intersect
    // with the edge of the triangle and populate the appropriate information.
    // Return true if you don't need to test anything else, false otherwise
    // Only test if triangle edge is not coplanar with the segment.  If it is
    // coplanar, it follows that one of the other edges of the triangle will
    // find the intersection (i.e. call exactIntersection with "false")
    if (exactIntersection(s0,s1,P,tA,tB,false,IntPt))
      {
	if (!IntPt)return true;
	Pts.push_back(*IntPt);
	if (Pts.size() == 2)
	  {
	    *IntPt = (Pts[0] + Pts[1]) * 0.5;
	    return true;
	  }
      }

    return false;
  }

  static bool testEdgeAlmostCoplanar(Vector3d const &s0,
                                     Vector3d const &s1,
                                     Vector3d const &P0,
                                     Vector3d const &P1,
                                     Vector3d const &tA,
                                     Vector3d const &tB,
                                     Vector3d &IntPt,
                                     std::vector< Vector3d > &Pts)
  {
    // alexr, 1-24-2014
    // This is essentially the same as testEdge, execpt it forms two triangles
    // with two different out of plane points and then tries to intersect each
    // of them, keeping at most one intersection point.
    //
    // IMPORTANT: this routine returns adds at most one point to the Pts vector...
    //   if the segment goes right through the edge, we don't want to count that
    //   intersection twice!
    //
    if (exactIntersection(s0,s1,P0,tA,tB,false,&IntPt))
      {
        Pts.push_back(IntPt);
        if (Pts.size() == 2)
          {
            IntPt = (Pts[0] + Pts[1]) * 0.5;
            return true;
          }
      }
    else if (exactIntersection(s0,s1,P1,tA,tB,false,&IntPt))
      {
         Pts.push_back(IntPt);
         if (Pts.size() == 2)
           {
             IntPt = (Pts[0] + Pts[1]) * 0.5;
             return true;
           }
      }
    return false;
  }


  static bool wrongSide(Vector3d const &t0,
			Vector3d const &t1,
			Vector3d const &t2,
			Vector3d const &s0,
			Vector3d const &s1)
  {
    double sign0 =
      orient3d((double*) &t0, (double*) &t1, (double*) &t2, (double*) &s0);
    if (sign0 >= 0)return false;
    double sign1 =
      orient3d((double*) &t0, (double*) &t1, (double*) &t2, (double*) &s1);
    if (sign1 >= 0)return false;
    return true;
  }

  static bool coplanarIntersection(Vector3d const &t0,
				   Vector3d const &t1,
				   Vector3d const &t2,
				   Vector3d const &s0,
				   Vector3d const &s1,
				   Vector3d *IntPt,
				   bool infinitePlane)
  {
    /* If you get here, you have the triangle and the segment in one plane!
     * 1) Lets find an out of plane point P.  If you cannot find it, you must
     *    be a colinear triangle -- return false
     *    If you find it, make sure that P is on the positive side of triangle
     * 2) Given P, see if end point 1 of the segment is inside triangle.  The
     *    way to do it is to form a segment with end point 1 and P and test for
     *    int
     * 3) Do the same for end point 2
     * 4) If neither point is inside, the only way you can have intersection is
     *    if any of the edges of the triangle intersects with the triangle
     *    containing the segment and point P.  Be careful here because you can
     *    have the case where the segment is colinear with a triangle edge, in
     *    which case, you get into another coplanar intersection test, which
     *    can become an infinite recursion!
     * 5) If there is intersection, you must have found 2 points.  If you need
     *    the location, return the average.
     *
     * NOTE:  Although this routine calls "exactIntersection", there is no risk
     *        of infinite recursion since every time we call it, we know that
     *        the triangle and segment provided to it will always be non-planar
     */

    // If you are testing intersection with an infinite plane, the result is
    // always true because we already lie on it!  Return the midpoint...
    if (infinitePlane)
      {
	if (IntPt)
	  *IntPt = (s0 + s1) * 0.5;
	return true;
      }

    Vector3d P;
    if (!findOutOfPlanePoint(t0,t1,t2,s0,s1,P))
      return false; // colinear triangle -- don't bother with intersection

    // The idea is that more often than not, we will not be intersecting.  So
    // lets do some work to filter out the non-intersecting guys first.  If we
    // cannot filter them out, we will do more work.  What I am saying is that
    // we can comment out the following 3 lines without changing the result!
    if (wrongSide(t1,t0,P,s0,s1))return false;
    if (wrongSide(t2,t1,P,s0,s1))return false;
    if (wrongSide(t0,t2,P,s0,s1))return false;

    // Store intersection points if needed
    static std::vector< Vector3d > Pts;
    Pts.clear();

    // Check if end point s0 is inside the triangle
    if (exactIntersection(t0,t1,t2,s0,P,false))
      {
	if (!IntPt)return true;
	Pts.push_back(s0);
      }

    // Check if end point s1 is inside the triangle
    if (exactIntersection(t0,t1,t2,s1,P,false))
      {
	if (!IntPt)return true;
	Pts.push_back(s1);
      }

    // If both end points are inside, we are done
    if (IntPt && (Pts.size() == 2))
      {
	*IntPt = (Pts[0] + Pts[1]) * 0.5;
	return true;
      }

    // Test each triangle edge with triangle formed using s0,s1,P
    if (testEdge(s0,s1,P,t0,t1,IntPt,Pts))return true;
    if (testEdge(s0,s1,P,t1,t2,IntPt,Pts))return true;
    if (testEdge(s0,s1,P,t2,t0,IntPt,Pts))return true;
    
    // If you get down here, you don't intersect.
    // Should we do a test to make sure we have exactly 2 ints?
    //    if(Pts.size() == 1)
    //      std::cout << "ONLY ONE INTERSECTION!!!" << std::endl;
    return false;
  }

  static bool isAlmostCoplanarIntersection(double volume,
                                           Vector3d const &t0,
                                           Vector3d const &t1,
                                           Vector3d const &t2)


  {
    double area = triangleArea(t0,t1,t2).mag();
    if (area < ALMOST_ZERO)
      return false;

    if (volume/area < 1.0e-12)
      return true;

    return false;
  }
  static void almostCoplanarIntersection(Vector3d const &t0,
                                         Vector3d const &t1,
                                         Vector3d const &t2,
                                         Vector3d const &s0,
                                         Vector3d const &s1,
                                         Vector3d &IntPt)
  {
    /*
     * alexr, 1-24-2014
     * This routine is uses an approach similar to the coplanarIntersection
     *   routine for finding the intersection point between a segment and
     *   a triangle.
     * This is being written to handle cases where the segment and triangle
     *   are nearly coplanar and the computation of the intersection point
     *   is inaccurate. The intention here is to prevent the intersection
     *   point computed from being far away from the triangle.
     *
     * The algorithm involves forming tetrahedra on either side of the triangle
     *   and itersecting the segment with the six exterior faces of this new
     *   polyhedron. We should find that either:
     *    * the two segment endpoints are inside the polyhedra,
     *    * two of the faces of these tetrahedra intersect the segment, or
     *    * one end point is inside and the segment intersects in one place.
     * In all cases, we find should be able to find two valid intersection points
     *  and take the midpoint of the two as the final result.
     * Rather than use the midpoint between these two, we will just return the
     *  first valid point, reducing the number of recursive calls we need to make.
     *
     */

    Vector3d P0;
    Vector3d P1;
    if (!findOutOfPlanePointsUsingNormal(t0,t1,t2,P0,P1))
      {
        if (!findOutOfPlanePoints(t0,t1,t2,s0,s1,P0,P1))
          {
            // if these methods all failed, I expect that the triangle is 
            // colinear which shouldn't be possibel if I entered this 
            // function... return the segment midpoint so something
            // is in IntPt            
             IntPt = (s0 + s1) * 0.5;
            return; 
          }
      }

    // Store intersection points if needed
    static std::vector< Vector3d > Pts;
    Pts.clear();

    if (exactIntersection(t0,t1,t2,s0,P0,false) ||
        exactIntersection(t0,t1,t2,s0,P1,false))
      {
        Pts.push_back(s0);
      }

    if (exactIntersection(t0,t1,t2,s1,P0,false) ||
        exactIntersection(t0,t1,t2,s1,P1,false))
      {
        Pts.push_back(s1);
      }

    // If both end points are inside, we are done
    if (Pts.size() == 2)
      {
        IntPt = (Pts[0] + Pts[1]) * 0.5;
        return;
      }

    // Test each triangle edge with triangle formed using s0,s1,P
    if (testEdgeAlmostCoplanar(s0,s1,P0,P1,t0,t1,IntPt,Pts))
      return;
    if (testEdgeAlmostCoplanar(s0,s1,P0,P1,t1,t2,IntPt,Pts))
      return;
    if (testEdgeAlmostCoplanar(s0,s1,P0,P1,t2,t0,IntPt,Pts))
      return;
  }

  // alexr, 1-24-2014, adding the almostCoplanarIntersection
  //                   case below to improve results in cases
  //                   where the near coplanarity was causing
  //                   the intersection point to inaccurate
  static bool exactIntersection(Vector3d const &t0,
                                Vector3d const &t1,
                                Vector3d const &t2,
                                Vector3d const &s0,
                                Vector3d const &s1,
                                bool testCoplanar,
                                Vector3d *IntPt,
                                bool infinitePlane)
  {
    double *A = (double*) &t0(0);
    double *B = (double*) &t1(0);
    double *C = (double*) &t2(0);
    double *D = (double*) &s0(0);
    double *E = (double*) &s1(0);

    double ABCE = orient3d(A,B,C,E);
    double ABCD = orient3d(A,B,C,D);

    if (ABCE < 0.0 || ABCD > 0.0) {
        double* tmpp;
        double  tmp;
        tmpp = E; E = D; D = tmpp;
        tmp = ABCE; ABCE = ABCD; ABCD = tmp;
    }
    if (ABCE < 0.0 || ABCD > 0.0) return false;

    // Treat special case of coplanar intersection
    if (ABCD == ABCE) // the only way they are equal is if they are both 0!
      {
        if (!testCoplanar)return false;
        return coplanarIntersection(t0,t1,t2,s0,s1,IntPt,infinitePlane);
      }

    if (!infinitePlane) // only if you are not testing int with infinite plane
      {
        if ((orient3d (A,D,C,E)) < 0.0) return false;
        if ((orient3d (A,B,D,E)) < 0.0) return false;
        if ((orient3d (B,C,D,E)) < 0.0) return false;
      }

    // At this point, we know we have non-coplanar intersection
    if (IntPt)
      {
        // if !testCoplanar, avoid the almostCoplanarIntersection cases
        //  to guarantee that we do recurse infinitely
        double volume = ABCE - ABCD;
        if (testCoplanar && volume < 1.0e-20 &&
            isAlmostCoplanarIntersection(volume,t0,t1,t2))
           {
             almostCoplanarIntersection(t0,t1,t2,s0,s1,*IntPt);
           }
        else
          {
            double c = ABCE/volume;
            (*IntPt)(0) = E[0] + c*(D[0] - E[0]);
            (*IntPt)(1) = E[1] + c*(D[1] - E[1]);
            (*IntPt)(2) = E[2] + c*(D[2] - E[2]);
          }

//        Triangle3 triangle(t0,t1,t2);
//        double distSq = triangle.getDistanceSquared(*IntPt);
//        if (distSq > 1.0e-20)
//          {
//            std::cout << "Bad Intersection Point, distSq=" << distSq
//                << ", ABCE=" << ABCE << ", ABCD=" << ABCD
//                << ", testCoplanar=" << testCoplanar << std::endl;
//          }
      }
    return true;
  }


  bool Triangle3::testExactIntersection(Segment3 const &segment,
					bool infinitePlane) const
  {
    return exactIntersection(this->getCoordinate(0),
			     this->getCoordinate(1),
			     this->getCoordinate(2),
			     segment.getCoordinate(0),
			     segment.getCoordinate(1),
			     true,
			     0,
			     infinitePlane);
  }


  bool Triangle3::findExactIntersection(Segment3 const &segment,
					Vector3d &IntPt,
					bool infinitePlane) const
  {
    return exactIntersection(this->getCoordinate(0),
			     this->getCoordinate(1),
			     this->getCoordinate(2),
			     segment.getCoordinate(0),
			     segment.getCoordinate(1),
			     true,
			     &IntPt,
			     infinitePlane);
  }




  // Exact intersection routines (with sos)
  // alexr, 8-2013, segmentOrientation parameter tells
  //                which way the segment is oriented with
  //                respect to the triangle normals...
  //                if the surface if oriented outward, this
  //                tells us if D or E is inside the surface
  // alexr, 1-24-2014, adding the almostCoplanarIntersection
  //                   case below to improve results in cases
  //                   where the near coplanarity was causing
  //                   the intersection point to inaccurate
  static bool exactIntersection_sos(Point3 const &A,
                                    Point3 const &B,
                                    Point3 const &C,
                                    Point3 const &D,
                                    Point3 const &E,
                                    Vector3d *IntPt=0,
                                    bool infinitePlane=false,
                                    bool * segmentOrientation=0)
  {
    Point3 const *Eptr = &E;
    Point3 const *Dptr = &D;
    double ABCEval;
    double ABCDval;
    double ABCE = orient3d_sos (A, B, C, *Eptr, &ABCEval);
    double ABCD = orient3d_sos (A, B, C, *Dptr, &ABCDval);

    if (ABCE < 0.0 || ABCD > 0.0) {
        Point3 const *tmpp = Eptr;  Eptr = Dptr;        Dptr = tmpp;
        double tmp = ABCE;          ABCE = ABCD;        ABCD = tmp;
        tmp = ABCEval;              ABCEval = ABCDval;  ABCDval = tmp;
        if (segmentOrientation)
          *segmentOrientation = false;
    } else if (segmentOrientation) {
        *segmentOrientation = true;
    }
    if (ABCE < 0.0 || ABCD > 0.0) return false;

    if (!infinitePlane) // only if you are not testing int with infinite plane
      {
        if ((orient3d_sos (A, *Dptr, C, *Eptr)) < 0.0) return false;
        if ((orient3d_sos (A, B, *Dptr, *Eptr)) < 0.0) return false;
        if ((orient3d_sos (B, C, *Dptr, *Eptr)) < 0.0) return false;
      }

    if (!IntPt)
      return true;

    // Now we need to find the intersection point
    if (ABCDval != ABCEval)
      {
        double volume = ABCEval - ABCDval;
        if (volume < 1.0e-20 &&
            isAlmostCoplanarIntersection(volume,A.position(),B.position(),C.position()))
          {
            almostCoplanarIntersection(A.position(),B.position(),C.position(),
                                       D.position(),E.position(),*IntPt);
          }
        else
          {
            double c = ABCEval / volume;
            *IntPt = Eptr->position() + c*(Dptr->position() - Eptr->position());
          }
      }
    else // coplanar -- need to do more work to find exact point of int!
      {
        // I don't know if this can ever return false, but if it does, just
        // return midpoint of segment as the intersection location
        if (!coplanarIntersection(A.position(),B.position(),C.position(),
                                  D.position(),E.position(),IntPt,
                                  infinitePlane))
          *IntPt = 0.5*(D.position() + E.position());
      }

//    Triangle3 triangle(A.position(), B.position(), C.position());
//    double distSq = triangle.getDistanceSquared(*IntPt);
//    if (distSq > 1.0e-20)
//      {
//        std::cout << "Bad Intersection Point, distSq=" << distSq
//            << ", ABCEval=" << ABCEval << ", ABCDval=" << ABCDval << std::endl;
//      }

    return true;
  }


  bool Triangle3::testExactIntersection_sos(Segment3 const &segment,
					    bool infinitePlane) const
  {
    return exactIntersection_sos(this->getPoint(0),
				 this->getPoint(1),
				 this->getPoint(2),
				 segment.getPoint(0),
				 segment.getPoint(1),
				 0,
				 infinitePlane);
  }




  bool Triangle3::findExactIntersection_sos(Segment3 const &segment,
					    Vector3d &IntPt,
					    bool infinitePlane,
					    bool * segmentOrientation) const
  {
    return exactIntersection_sos(this->getPoint(0),
				 this->getPoint(1),
				 this->getPoint(2),
				 segment.getPoint(0),
				 segment.getPoint(1),
				 &IntPt,
				 infinitePlane,
				 segmentOrientation);
  }


#if 0 // old altorithm (had tolerance issues)...
  // use the Moller/Trumbore algorithm from the paper:
  //    "Fast, Minimum Storage Ray-Triangle Intersection"
  // following the non-culling branch
  bool Triangle3::testIntersection(Ray3 const &ray,
				   double *distance) const
  {
    // fixme -- should this be adjustable?
#if 1
    static double DETERMINANT_TOL = 1.0e-9;
#else
    static double DETERMINANT_TOL = ALMOST_ZERO;
#endif

    /* find vectors for two edges sharing vertex 1 */
    Vector<3, double> edge1 = this->getPoint(1).position();
    edge1 -= this->getPoint(0).position();
    Vector<3, double> edge2 = this->getPoint(2).position();
    edge2 -= this->getPoint(0).position();

    Vector3d const &dir = ray.getDirection();
        
    /* begin calculating determinant - also used to calculate U parameter*/
    Vector3d pvec = dir.cross(edge2);

    /* if determinant is near zero, ray lies in plane of triangle */
    double det = edge1.dot(pvec);
    if (det > -DETERMINANT_TOL && det < DETERMINANT_TOL)
      return false;
    double inv_det = 1.0 / det;

    /* calculate distance from vertex 1 to ray origin */
    Vector3d tvec(ray.getOrigin() - this->getCoordinate(0));

    /* calculate U parameter and test bounds */
    double u = tvec.dot(pvec) * inv_det;
    if (u < 0.0 || u > 1.0)
      return false;

    /* prepare to test V parameter */
    Vector3d qvec = tvec.cross(edge1);

    /* calculate V parameter and test bounds */
    double v = dir.dot(qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0)
      return false;

    /*! line segment passing through origin intersects triangle, but need
      to ensure that ray intersects triangle in positive direction.  Compute
      distance to check this.
    */

    double t = edge2.dot(qvec) * inv_det;
    if (t < 0.0)
      return false;

    /* ray intersects triangle. */
    if (distance)
      *distance = t;

    return true;
  }
#else
  // The above approach caused problems when the triangle was very small and
  // falsely reported no intersection since we are comparing a dimensional
  // value with 1.e-9 (i.e. the triangle had not been normalized).  I am
  // cooking up a new way to test for Triangle/Ray intersection below
  bool Triangle3::testIntersection(Ray3 const &ray,
				   double *distance) const
  {
    /* I had a case that failed on linux 32 bit platform only (ccmp2575).  The
     * correct fix would have been to convert this to exact intersection once
     * and for all.  However, I saw that this slowed down the code by 30-40%
     * and did not have the time to find ways to improve it.  Also, the Ray3
     * class is really not written well (for example, we don't know if the
     * vector is normalized) and would need to be cleaned up for a formal
     * commit.  I found an alternate solution to the bug at hand (see below).
     * If you ever have to revisit this routine, it may be time to bite the
     * bullet.  The #if 0 code below is the exact precision version I tried...
     * ALY: 04/17/2008
     */
#if 0
    static BoundingBox bbox;
    bbox = this->getBoundingBox();
    bbox.expand(ray.getOrigin());
    double len = 100*bbox.longestSideLength();
    Vector3d p2(ray.getOrigin());
    p2 += len*ray.getDirection();
    Vector3d intPt;
    if (exactIntersection(this->getCoordinate(0),
                         this->getCoordinate(1),
                         this->getCoordinate(2),
                         ray.getOrigin(),
                         p2,false,&intPt))
    {
      if (distance)
	*distance = intPt.distance(ray.getOrigin());
      return true;
    }
    return false;
#endif
    // Triangle can only intersect ray if ray origin is above triangle and
    // pointed down, or if origin is below triangle and pointed up.  Lets test
    // that first
    double dot = this->getNormal().dot(ray.getDirection());
    double sign = 1.0;

    // Assume no intersection if ray parallel to triangle plane
    if (fabs(dot) <= ALMOST_ZERO)return false;

    Vector3d vecToOrigin(ray.getOrigin() - this->getCoordinate(0));
    double dot2 = this->getNormal().dot(vecToOrigin);
    if (dot > 0.0) // origin must be below triangle
      {
	if (dot2 > 0.0)return false;
      }
    else // origin must be above triangle
      {
	if (dot2 < 0.0)return false;
	sign = -1.0;
      }

    /* Found a case on linux32 where a ray was parallel and coplanar but dot
     * above was 1.e-20.  We were finding intersections when we shouldn't have.
     * As a fix, I am checking to see if dot2 is almost zero (ray origin is
     * coplanar with triangle).  If so, then instead of using the ray normal
     * to compute the parallel piped below, I will use the triangle normal
     * because all we really need to check in this special case is if the
     * origin falls inside the triangle.
     * ALY: 04/17/2008
     * Note: I am not too thrilled about this solution as the vectors use to
     * compute dot2 have not necessarily been normalized.  However, since we
     * are comparing with a VERY small number, chances are this is okay!
     * ALY: 04/18/2008
     */
    Vector3d prismDir(ray.getDirection());
    if (fabs(dot2) <= ALMOST_ZERO)
      prismDir = this->getNormal();

    // Now you know stuff is on the right half plane.  We simply need to see if
    // it falls inside the triangle.  Do that by creating 3 plane vectors (one
    // for each edge) and making sure that the point is on the inside of all of
    // them.
    Vector3d cross = prismDir.cross(getEdge(0).getDirection());
    if (sign*cross.dot(vecToOrigin) >= 0.0)
      {
	cross = prismDir.cross(getEdge(1).getDirection());
	Vector3d vec(ray.getOrigin() - getCoordinate(1));
	if (sign*cross.dot(vec) >= 0.0)
	  {
	    cross = prismDir.cross(getEdge(2).getDirection());
	    vec   = ray.getOrigin() - getCoordinate(2);
	    if (sign*cross.dot(vec) >= 0.0) // ray intersects triangle
	      {
		if (distance)
		  {
		    double ratio = -dot2 / dot;
		    if (ratio < 0.0)ratio = 0.0; // should not happen
		    *distance = ray.getDirection().mag()*ratio;
		  }
		return true;
	      }
	  }
      }

    return false;
  }
#endif




/*
What follows is intersection tests for line-triangle. It's derived
from the above ray-intersection code, with sections that never run removed.
*/

  bool Triangle3::testIntersection(Line3 const &line,
				   double *distance) const // signed distance!
  {
    //Triangle and line do not intersect if the line is parallel to the plane
    double dot = this->getNormal().dot(line.getDirection());
    double sign = 1.0;
    if (fabs(dot) <= ALMOST_ZERO)return false;

    Vector3d vecToOrigin(line.getOrigin() - this->getCoordinate(0));
    double dot2 = this->getNormal().dot(vecToOrigin);

    if (dot < 0.0) 
      {
	sign = -1.0;
      }

    /* Found a case on linux32 where a ray was parallel and coplanar but dot
     * above was 1.e-20.  We were finding intersections when we shouldn't have.
     * As a fix, I am checking to see if dot2 is almost zero (ray origin is
     * coplanar with triangle).  If so, then instead of using the ray normal
     * to compute the parallel piped below, I will use the triangle normal
     * because all we really need to check in this special case is if the
     * origin falls inside the triangle.
     * ALY: 04/17/2008
     * Note: I am not too thrilled about this solution as the vectors use to
     * compute dot2 have not necessarily been normalized.  However, since we
     * are comparing with a VERY small number, chances are this is okay!
     * ALY: 04/18/2008
     */
    Vector3d prismDir(line.getDirection());
    if (fabs(dot2) <= ALMOST_ZERO)
      prismDir = this->getNormal();
    
    // Now you know stuff is on the right half plane.  We simply need to see if
    // it falls inside the triangle.  Do that by creating 3 plane vectors (one
    // for each edge) and making sure that the point is on the inside of all of
    // them.
    Vector3d cross = prismDir.cross(getEdge(0).getDirection());
    if (sign*cross.dot(vecToOrigin) >= 0.0)
      {
	cross = prismDir.cross(getEdge(1).getDirection());
	Vector3d vec(line.getOrigin() - getCoordinate(1));
	if (sign*cross.dot(vec) >= 0.0)
	  {
	    cross = prismDir.cross(getEdge(2).getDirection());
	    vec   = line.getOrigin() - getCoordinate(2);
	    if (sign*cross.dot(vec) >= 0.0) // ray intersects triangle
	      {
		if (distance)
		  {
		    double ratio = -dot2 / dot;
		    /* This was badly copied from the triangle-ray
		     * intersection!  In the case of triangle-line intersection
		     * the ratio can definitely be negative because that just
		     * implies that the intersection is in the other direction
		     * with respect to the origin and the line direction!
		     *
		     * NOTE:  This means that we are returning a signed
		     *        distance!  User-beware!
		     *
		     * This is the line that shouldn't exist!
		     * if(ratio < 0.0)ratio = 0.0; // should not happen
		     *
		     * ALY: 04/13/2009
		     */
		    *distance = line.getDirection().mag()*ratio;
		  }
		return true;
	      }
	  }
      }

    return false;
  }


  bool
  Triangle3::
  testExactIntersection(Ray3 const &ray, double *distance, double *Dot) const
  {
    /*
     * This is essentially a refactoring of the non-exact coding in
     *
     *   bool Triangle3::testIntersection(Ray3 const &, double *) const
     *
     * using the identity
     *
     *   a.(b x c) \equiv orient3d(a,b,c,0)
     */
    Vector3d const &d = ray.getDirection();
    Vector3d const &e0 = getEdge(0).getDirection();
    Vector3d const &e1 = getEdge(1).getDirection();
    Vector3d const zero(0,0,0);

    double const dot = orient3d((double*)&d(0),
                                (double*)&e0(0),
                                (double*)&e1(0),
                                (double*)&zero(0));

    if (Dot)
      *Dot = 0.5*dot;

    // Assume no intersection if ray parallel to plane of triangle
    if (dot == 0)
      return false;

    Vector3d const &p = ray.getOrigin();
    Vector3d s(p - getCoordinate(0));
    double const dot2 = orient3d((double*)&s(0),
                                 (double*)&e0(0),
                                 (double*)&e1(0),
                                 (double*)&zero(0));

    // Ray is moving away from plane of triangle
    if ((dot > 0 && dot2 > 0) || (dot < 0 && dot2 < 0))
      return false;

    // Now you know stuff is on the right half plane.  We simply need to see if
    // it falls inside the triangle.  Do that by creating 3 plane vectors (one
    // for each edge) and making sure that the point is on the inside of all of
    // them.
    int const sign = dot > 0 ? 1 : -1;

    if (sign*orient3d((double*)&s(0),
                      (double*)&d(0),
                      (double*)&e0(0),
                      (double*)&zero(0)) < 0) return false;

    s = p - getCoordinate(1);
    if (sign*orient3d((double*)&s(0),
                      (double*)&d(0),
                      (double*)&e1(0),
                      (double*)&zero(0)) < 0) return false;

    s = p - getCoordinate(2);
    Vector3d const &e2 = getEdge(2).getDirection();
    if (sign*orient3d((double*)&s(0),
                      (double*)&d(0),
                      (double*)&e2(0),
                      (double*)&zero(0)) < 0) return false;

    // There is definitely an intersection
    if (distance)
      *distance = -dot2/dot;

    return true;
  }


  bool
  Triangle3::
  testExactIntersection(Line3 const &line, double *distance, double *Dot) const
  {
    /*
     * The same as
     *
     *   bool Triangle3::testExactIntersection(Ray3 const &, ...) const
     *
     * but with allowance for signed distance
     */
    Vector3d const &d = line.getDirection();
    Vector3d const &e0 = getEdge(0).getDirection();
    Vector3d const &e1 = getEdge(1).getDirection();
    Vector3d const zero(0,0,0);

    double const dot = orient3d((double*)&d(0),
                                (double*)&e0(0),
                                (double*)&e1(0),
                                (double*)&zero(0));

    if (Dot)
      *Dot = 0.5*dot;

    // Assume no intersection if line parallel to plane of triangle
    if (dot == 0)
      return false;

    Vector3d const &p = line.getOrigin();
    Vector3d s(p - getCoordinate(0));
    double const dot2 = orient3d((double*)&s(0),
                                 (double*)&e0(0),
                                 (double*)&e1(0),
                                 (double*)&zero(0));

    // Now you know stuff is on the right half plane.  We simply need to see if
    // it falls inside the triangle.  Do that by creating 3 plane vectors (one
    // for each edge) and making sure that the point is on the inside of all of
    // them.
    int const sign = dot > 0 ? 1 : -1;

    if (sign*orient3d((double*)&s(0),
                      (double*)&d(0),
                      (double*)&e0(0),
                      (double*)&zero(0)) < 0) return false;

    s = p - getCoordinate(1);
    if (sign*orient3d((double*)&s(0),
                      (double*)&d(0),
                      (double*)&e1(0),
                      (double*)&zero(0)) < 0) return false;

    s = p - getCoordinate(2);
    Vector3d const &e2 = getEdge(2).getDirection();
    if (sign*orient3d((double*)&s(0),
                      (double*)&d(0),
                      (double*)&e2(0),
                      (double*)&zero(0)) < 0) return false;

    // There is definitely an intersection
    if (distance)
      *distance = -dot2/dot;

    return true;
  }
}
