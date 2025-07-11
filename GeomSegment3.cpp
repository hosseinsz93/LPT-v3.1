#include "GeomSegment3.h"
#include "GeomBox3.h"
#include "GeomUtility.h"

namespace Geom
{
  void Segment3::computeUnitDirection() const
  {
    _unitdirection = getDirection();
    _length = _unitdirection.normalize();
    _hasunitdirection = true;
  }

  void Segment3::computeBoundingBox() const
  {
    _boundingBox.reset(_vertices[0].position());
    _boundingBox.expand(_vertices[1].position());
    _hasBoundingBox = true;
  }
  
  BoundingBox Segment3::getBoundingBox() const
  {
    if(!_hasBoundingBox)
      {
	computeBoundingBox();
      }
    return _boundingBox;
  }
  
  BoundingBox const& Segment3::getConstBoundingBox() const
  {
    if(!_hasBoundingBox)
      {
	computeBoundingBox();
      }
    return _boundingBox;
  }


  bool Segment3::intersectsBox(BoundingBox const &bbox) const
  {
    return Box3(bbox.min(), bbox.max()).testIntersection(*this);
  }
  

  bool Segment3::containsPoint(Vector3d const &point) const
  {
    return (getDistanceSquared(point) <= ALMOST_ZERO);
  }

  double Segment3::getDistanceSquared (Segment3 const &rkSeg1, 
				       double* pfSegP0, double* pfSegP1) const
  {
    Vector3d const zero(0.0,0.0,0.0);
    Vector3d kDiff(this->getCoordinate(0) - rkSeg1.getCoordinate(0));
    double fA00 = this->getDirection().mag2();
    
    double fA01 = Vector3d(zero-(this->getDirection())).dot(rkSeg1.getDirection());
    double fA11 = rkSeg1.getDirection().mag2();
    double fB0 = kDiff.dot(this->getDirection());
    double fC = kDiff.mag2();
    double fDet = fabs(fA00*fA11-fA01*fA01);
    double fB1, fS, fT, fSqrDist, fTmp;

    static double DETERMINANT_TOL = 1.0e-9;//should it be adjustable?

    if ( fDet >= DETERMINANT_TOL )
      {
        // line segments are not parallel
        fB1 = Vector3d(zero-kDiff).dot(rkSeg1.getDirection());
        fS = fA01*fB1-fA11*fB0;
        fT = fA01*fB0-fA00*fB1;
        if ( fS >= (double)0.0 )
          {
            if ( fS <= fDet )
              {
                if ( fT >= (double)0.0 )
                  {
                    if ( fT <= fDet )  // region 0 (interior)
                      {
                        // minimum at two interior points of 3D lines
                        double fInvDet = ((double)1.0)/fDet;
                        fS *= fInvDet;
                        fT *= fInvDet;
                        fSqrDist = fS*(fA00*fS+fA01*fT+((double)2.0)*fB0) +
                          fT*(fA01*fS+fA11*fT+((double)2.0)*fB1)+fC;
                      }
                    else  // region 3 (side)
                      {
                        fT = (double)1.0;
                        fTmp = fA01+fB0;
                        if ( fTmp >= (double)0.0 )
                          {
                            fS = (double)0.0;
                            fSqrDist = fA11+((double)2.0)*fB1+fC;
                          }
                        else if ( -fTmp >= fA00 )
                          {
                            fS = (double)1.0;
                            fSqrDist = fA00+fA11+fC+((double)2.0)*(fB1+fTmp);
                          }
                        else
                          {
                            fS = -fTmp/fA00;
                            fSqrDist = fTmp*fS+fA11+((double)2.0)*fB1+fC;
                          }
                      }
                  }
                else  // region 7 (side)
                  {
                    fT = (double)0.0;
                    if ( fB0 >= (double)0.0 )
                      {
                        fS = (double)0.0;
                        fSqrDist = fC;
                      }
                    else if ( -fB0 >= fA00 )
                      {
                        fS = (double)1.0;
                        fSqrDist = fA00+((double)2.0)*fB0+fC;
                      }
                    else
                      {
                        fS = -fB0/fA00;
                        fSqrDist = fB0*fS+fC;
                      }
                  }
              }
            else
              {
                if ( fT >= (double)0.0 )
                  {
                    if ( fT <= fDet )  // region 1 (side)
                      {
                        fS = (double)1.0;
                        fTmp = fA01+fB1;
                        if ( fTmp >= (double)0.0 )
                          {
                            fT = (double)0.0;
                            fSqrDist = fA00+((double)2.0)*fB0+fC;
                          }
                        else if ( -fTmp >= fA11 )
                          {
                            fT = (double)1.0;
                            fSqrDist = fA00+fA11+fC+((double)2.0)*(fB0+fTmp);
                          }
                        else
                          {
                            fT = -fTmp/fA11;
                            fSqrDist = fTmp*fT+fA00+((double)2.0)*fB0+fC;
                          }
                      }
                    else  // region 2 (corner)
                      {
                        fTmp = fA01+fB0;
                        if ( -fTmp <= fA00 )
                          {
                            fT = (double)1.0;
                            if ( fTmp >= (double)0.0 )
                              {
                                fS = (double)0.0;
                                fSqrDist = fA11+((double)2.0)*fB1+fC;
                              }
                            else
                              {
                                fS = -fTmp/fA00;
                                fSqrDist = fTmp*fS+fA11+((double)2.0)*fB1+fC;
                              }
                          }
                        else
                          {
                            fS = (double)1.0;
                            fTmp = fA01+fB1;
                            if ( fTmp >= (double)0.0 )
                              {
                                fT = (double)0.0;
                                fSqrDist = fA00+((double)2.0)*fB0+fC;
                              }
                            else if ( -fTmp >= fA11 )
                              {
                                fT = (double)1.0;
                                fSqrDist = fA00+fA11+fC+
                                  ((double)2.0)*(fB0+fTmp);
                              }
                            else
                              {
                                fT = -fTmp/fA11;
                                fSqrDist = fTmp*fT+fA00+((double)2.0)*fB0+fC;
                              }
                          }
                      }
                  }
                else  // region 8 (corner)
                  {
                    if ( -fB0 < fA00 )
                      {
                        fT = (double)0.0;
                        if ( fB0 >= (double)0.0 )
                          {
                            fS = (double)0.0;
                            fSqrDist = fC;
                          }
                        else
                          {
                            fS = -fB0/fA00;
                            fSqrDist = fB0*fS+fC;
                          }
                      }
                    else
                      {
                        fS = (double)1.0;
                        fTmp = fA01+fB1;
                        if ( fTmp >= (double)0.0 )
                          {
                            fT = (double)0.0;
                            fSqrDist = fA00+((double)2.0)*fB0+fC;
                          }
                        else if ( -fTmp >= fA11 )
                          {
                            fT = (double)1.0;
                            fSqrDist = fA00+fA11+fC+((double)2.0)*(fB0+fTmp);
                          }
                        else
                          {
                            fT = -fTmp/fA11;
                            fSqrDist = fTmp*fT+fA00+((double)2.0)*fB0+fC;
                          }
                      }
                  }
              }
          }
        else 
          {
            if ( fT >= (double)0.0 )
              {
                if ( fT <= fDet )  // region 5 (side)
                  {
                    fS = (double)0.0;
                    if ( fB1 >= (double)0.0 )
                      {
                        fT = (double)0.0;
                        fSqrDist = fC;
                      }
                    else if ( -fB1 >= fA11 )
                      {
                        fT = (double)1.0;
                        fSqrDist = fA11+((double)2.0)*fB1+fC;
                      }
                    else
                      {
                        fT = -fB1/fA11;
                        fSqrDist = fB1*fT+fC;
                      }
                  }
                else  // region 4 (corner)
                  {
                    fTmp = fA01+fB0;
                    if ( fTmp < (double)0.0 )
                      {
                        fT = (double)1.0;
                        if ( -fTmp >= fA00 )
                          {
                            fS = (double)1.0;
                            fSqrDist = fA00+fA11+fC+((double)2.0)*(fB1+fTmp);
                          }
                        else
                          {
                            fS = -fTmp/fA00;
                            fSqrDist = fTmp*fS+fA11+((double)2.0)*fB1+fC;
                          }
                      }
                    else
                      {
                        fS = (double)0.0;
                        if ( fB1 >= (double)0.0 )
                          {
                            fT = (double)0.0;
                            fSqrDist = fC;
                          }
                        else if ( -fB1 >= fA11 )
                          {
                            fT = (double)1.0;
                            fSqrDist = fA11+((double)2.0)*fB1+fC;
                          }
                        else
                          {
                            fT = -fB1/fA11;
                            fSqrDist = fB1*fT+fC;
                          }
                      }
                  }
              }
            else   // region 6 (corner)
              {
                if ( fB0 < (double)0.0 )
                  {
                    fT = (double)0.0;
                    if ( -fB0 >= fA00 )
                      {
                        fS = (double)1.0;
                        fSqrDist = fA00+((double)2.0)*fB0+fC;
                      }
                    else
                      {
                        fS = -fB0/fA00;
                        fSqrDist = fB0*fS+fC;
                      }
                  }
                else
                  {
                    fS = (double)0.0;
                    if ( fB1 >= (double)0.0 )
                      {
                        fT = (double)0.0;
                        fSqrDist = fC;
                      }
                    else if ( -fB1 >= fA11 )
                      {
                        fT = (double)1.0;
                        fSqrDist = fA11+((double)2.0)*fB1+fC;
                      }
                    else
                      {
                        fT = -fB1/fA11;
                        fSqrDist = fB1*fT+fC;
                      }
                  }
              }
          }
      }
    else
      {
        // line segments are parallel
        if ( fA01 > (double)0.0 )
          {
            // direction vectors form an obtuse angle
            if ( fB0 >= (double)0.0 )
              {
                fS = (double)0.0;
                fT = (double)0.0;
                fSqrDist = fC;
              }
            else if ( -fB0 <= fA00 )
              {
                fS = -fB0/fA00;
                fT = (double)0.0;
                fSqrDist = fB0*fS+fC;
              }
            else
              {
                fB1 = Vector3d(zero-kDiff).dot(rkSeg1.getDirection());
                fS = (double)1.0;
                fTmp = fA00+fB0;
                if ( -fTmp >= fA01 )
                  {
                    fT = (double)1.0;
                    fSqrDist = fA00+fA11+fC+((double)2.0)*(fA01+fB0+fB1);
                  }
                else
                  {
                    fT = -fTmp/fA01;
                    fSqrDist = fA00+((double)2.0)*fB0+fC+fT*(fA11*fT+
                                                             ((double)2.0)*(fA01+fB1));
                  }
              }
          }
        else
          {
            // direction vectors form an acute angle
            if ( -fB0 >= fA00 )
              {
                fS = (double)1.0;
                fT = (double)0.0;
                fSqrDist = fA00+((double)2.0)*fB0+fC;
              }
            else if ( fB0 <= (double)0.0 )
              {
                fS = -fB0/fA00;
                fT = (double)0.0;
                fSqrDist = fB0*fS+fC;
              }
            else
              {
                fB1 = Vector3d(zero-kDiff).dot(rkSeg1.getDirection());
                fS = (double)0.0;
                if ( fB0 >= -fA01 )
                  {
                    fT = (double)1.0;
                    fSqrDist = fA11+((double)2.0)*fB1+fC;
                  }
                else
                  {
                    fT = -fB0/fA01;
                    fSqrDist = fC+fT*(((double)2.0)*fB1+fA11*fT);
                  }
              }
          }
      }

    if ( pfSegP0 )
      *pfSegP0 = fS;

    if ( pfSegP1 )
      *pfSegP1 = fT;

    return fabs(fSqrDist);
  }






  /* I (Aly) had precision problems with the above routine when two segments
   * penetrated each other (returned non-zero distance).  I have rewritten it
   * my way and it seems to work.  I will just call it something different and
   * if it seems to be okay for a while, we can just remove the above routine.
   * 3/25/2005
   */
  double Segment3::getDistanceSquared2(Segment3 const &segment2,
				       double* t1, double* t2) const
  {
    // Return distance from point 1 to segment 2 if segment 1 is 0 length
    double mag1sq = this->getDirection().dot(this->getDirection());
    if (mag1sq < ALMOST_ZERO)
      {
	if (t1) *t1 = 0.0;
	return segment2.getDistanceSquared(this->getCoordinate(0),0,t2);
      }

    // Return distance from point 2 to segment 1 if segment 2 is 0 length
    double mag2sq = segment2.getDirection().dot(segment2.getDirection());
    if (mag2sq < ALMOST_ZERO)
      {
	if (t2) *t2 = 0.0;
	return this->getDistanceSquared(segment2.getCoordinate(0),0,t1);
      }

    double dot = this->getDirection().dot(segment2.getDirection());

    double cosSq = (dot*dot)/(mag1sq*mag2sq);

    // Return minimum of either endpoint on segment2 to segment1 if parallel
    if ((1.0-cosSq) < 1.e-8) // parallel
      {
	if (mag1sq > mag2sq)
	  {
	    double d1 = this->getDistanceSquared(segment2.getCoordinate(0),0,t1);
	    double d2 = this->getDistanceSquared(segment2.getCoordinate(1),0,t1);
	    
	    if (d1 < d2)
	      {
		if (t2) *t2 = 0.0;
		return d1;
	      }
	    else
	      {
		if (t2) *t2 = 1.0;
		return d2;
	      }
	  }
	else
	  {
	    double d1 = segment2.getDistanceSquared(this->getCoordinate(0),0,t2);
	    double d2 = segment2.getDistanceSquared(this->getCoordinate(1),0,t2);
	    
	    if (d1 < d2)
	      {
		if (t1) *t1 = 0.0;
		return d1;
	      }
	    else
	      {
		if (t1) *t1 = 1.0;
		return d2;
	      }
	  }
      }

    // Segments are not parallel.  The shortest distance between the two is
    // orthogonal to both of them
    Vector3d orthog = this->getDirection().cross(segment2.getDirection());
    // Create a plane with segment 1 and the orthogonal and intersect segment2
    // with it...
    Vector3d plane = this->getDirection().cross(orthog);
    Vector3d vec1(segment2.getCoordinate(0) - this->getCoordinate(0));
    Vector3d vec2(segment2.getCoordinate(1) - this->getCoordinate(0));
    double dot1 = vec1.dot(plane);
    double dot2 = vec2.dot(plane);

    // seg 2 on +ve side of plane
    if (dot1 >= 0.0 && dot2 >= 0.0)
      {
	if (dot1 < dot2)
	  {
	    if (t2) *t2 = 0.0;
	    return this->getDistanceSquared(segment2.getCoordinate(0),0,t1);
	  }
	else
	  {
	    if (t2) *t2 = 1.0;
	    return this->getDistanceSquared(segment2.getCoordinate(1),0,t1);
	  }
      }

    // seg 2 on -ve side of plane
    if (dot1 <= 0.0 && dot2 <= 0.0)
      {
	if (dot1 < dot2)
	  {
	    if (t2) *t2 = 1.0;
	    return this->getDistanceSquared(segment2.getCoordinate(1),0,t1);
	  }
	else
	  {
	    if (t2) *t2 = 0.0;
	    return this->getDistanceSquared(segment2.getCoordinate(0),0,t1);
	  }
      }

    // seg 2 crosses plane
    double ratio = dot1/(dot1-dot2); // signs work themselves out...
    Vector3d P(segment2.getCoordinate(0) +
		       segment2.getDirection()*ratio);
    if (t2) *t2 = ratio;
    return this->getDistanceSquared(P,0,t1);
  }

  double Segment3::getDistanceSquared(Vector3d const &querypt,
				      Vector3d *closestpoint,
				      double *param) const
  {
    // distance from point to segment
    Vector3d diff(querypt - this->getCoordinate(0));
    Vector3d dir(this->getCoordinate(1) - this->getCoordinate(0));
    double t = diff.dot(dir);

    if (t <= 0.0) // point is behind first endpoint of edge
      {
	if (closestpoint)
	  *closestpoint = this->getCoordinate(0);
	if (param)
	  *param = 0.0;
	return diff.mag2();
      }

    double len_squared = dir.mag2();
    if (t >= len_squared) // point is infront of second endpoint of edge
      {
	if (closestpoint)
	  *closestpoint = this->getCoordinate(1);
	if (param)
	  *param = 1.0;
	diff -= dir;
	return diff.mag2();
      }

    // point projects between two edge endpoints
    t /= len_squared;
    diff -= dir*t;
    if (closestpoint)
      *closestpoint = this->getCoordinate(0) + dir*t;
    if (param)
      *param = t;
    return diff.mag2();
  }
}
