#include "GeomPlane3.h"

// #include "platform/OS.h"

#include "GeomLine3.h"
#include "GeomSegment3.h"
#include "GeomUtility.h"

namespace Geom
{
  void Plane3::computeEquation() const
  {
    _eqn[0] = _normal(0);
    _eqn[1] = _normal(1);
    _eqn[2] = _normal(2);
    _eqn[3] = _normal.dot(_point);
    _haseqn = true;
  }


  bool Plane3::findIntersection(Plane3 const &plane2,
				Line3 &intLine) const
  {
    // Taken from www.magic-software.com

    // If Cross(N0,N1) is zero, then either planes are parallel and separated
    // or the same plane.  In both cases, 'false' is returned.  Otherwise,
    // the intersection line is
    //
    //   L(t) = t*Cross(N0,N1) + c0*N0 + c1*N1
    //
    // for some coefficients c0 and c1 and for t any real number (the line
    // parameter).  Taking dot products with the normals,
    //
    //   d0 = Dot(N0,L) = c0*Dot(N0,N0) + c1*Dot(N0,N1)
    //   d1 = Dot(N1,L) = c0*Dot(N0,N1) + c1*Dot(N1,N1)
    //
    // which are two equations in two unknowns.  The solution is
    //
    //   c0 = (Dot(N1,N1)*d0 - Dot(N0,N1)*d1)/det
    //   c1 = (Dot(N0,N0)*d1 - Dot(N0,N1)*d0)/det
    //
    // where det = Dot(N0,N0)*Dot(N1,N1)-Dot(N0,N1)^2.

    // fixme -- should this be adjustable?
    static double DETERMINANT_TOL = 1.0e-9;

    double fN00 = this->getNormal().mag2();
    double fN01 = this->getNormal().dot(plane2.getNormal());
    double fN11 = plane2.getNormal().mag2();
    double fDet = fN00*fN11 - fN01*fN01;

    if (fabs(fDet) < DETERMINANT_TOL)
      return false;

    double fInvDet = ((double)1.0)/fDet;
    double fC0 = (fN11*this->getConstant() -
                  fN01*plane2.getConstant())*fInvDet;
    double fC1 = (fN00*plane2.getConstant() -
                  fN01*this->getConstant())*fInvDet;

    intLine.setDirection(this->getNormal().cross(plane2.getNormal()));
    intLine.setOrigin(Vector3d(this->getNormal()*fC0 + plane2.getNormal()*fC1));

    return true;
  }

  double Plane3::getSignedDistance(Vector3d const &point,
				   Vector3d *closestpoint) const
  {
    double const *eqn = this->getEquation();
    double distance = eqn[0]*point(0) + eqn[1]*point(1) + 
      eqn[2]*point(2) - eqn[3];
    if (closestpoint)
      {
        *closestpoint = point - this->getNormal()*distance;
      }
    return distance;
  }

  // this can probably be made faster 
  //returns -1 for segloc if the segment lies in the plane
  //otherwise segloc is the parametric distance from coordinate1 
  //that the intersection occurs i.e. v1 + (v2-v1)*segloc is the
  //point the intersection occurs where v1 and v2 are the vectors
  //indicating the endpoints of the segment
  bool Plane3::findIntersection(Segment3 const &seg,
				double* segloc) const
  {
    double dot1,dot2;
    Vector3d v1(seg.getCoordinate(0) - this->getPoint());
    Vector3d v2(seg.getCoordinate(1) - this->getPoint());
    
    dot1 = v1.dot(this->getNormal());
    dot2 = v2.dot(this->getNormal());
    
    if ((dot1 <= ALMOST_ZERO) && (dot1 >= -ALMOST_ZERO))
      {
	if ((dot2 <= ALMOST_ZERO) && (dot2 >= -ALMOST_ZERO))
	  {
	    *segloc = -1;//seg lies in the plane
	    return true;
	  }
	else
	  {
	      *segloc = 0;//endpoint v1 lies in plane
	    return true;
	  }
      }
    if ((dot2 <= ALMOST_ZERO) && (dot2 >= -ALMOST_ZERO))
      {
	if ((dot1 <= ALMOST_ZERO) && (dot1 >= -ALMOST_ZERO))
	  
	  {
	    *segloc = -1;//seg lies in the plane
	    return true;
	  }
	else
	  {
	      *segloc = 1.0;//endpoint v2 lies in plane
	      return true;
	  }
	
      }
    
    if (dot1*dot2 > 0.0)//both points on same side of plane so no intersection 
      {
	return false;
      }

    Vector3d v = seg.getDirection();
    //chrisgrevisit need to take square root - I can avoid this
    double denom = v.dot(this->getNormal());
    //in theory if we get here denom should never be 0 however, due to 
    //floating point error it can happen that dot1 and dot2 are very small
    //but larger than ALMOST_ZERO in such a way that denom comes out to be
    //0.  If this is true, we will count this as being in the plane.
    if(fabs(denom) < ALMOST_ZERO)
    {
      *segloc = -1;
      return true;
    }
    *segloc = -dot1/denom;
    return true;
  }

  bool Plane3::findIntersection(Line3 const &line,
				double* linloc) const
  {
    /* I didn't find this in the magic library so I am putting my own version.
     * At some point we really need to come back and make all of these
     * consistent with each other with regards to tolerances, etc.
     */

    // If line is parallel to plane, return false
    double denom = this->getNormal().dot(line.getDirection());
    if (denom > -ALMOST_ZERO && denom < ALMOST_ZERO)
      return false;

    Vector3d vec(this->getPoint() - line.getOrigin());
    double num = this->getNormal().dot(vec);

    *linloc = num/denom;

    return true;
  }
}
