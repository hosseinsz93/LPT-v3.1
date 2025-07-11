#include "GeomLine3.h"

#include "GeomPoint3.h"
#include "GeomSegment3.h"
#include "GeomUtility.h"

namespace Geom
{
  double Line3::getDistanceSquared(Vector3d const &point,
				   double* linloc) const
  {
    Vector3d kDiff(point - getOrigin());
    double fSqrLen = getDirection().mag2();
    double fT = kDiff.dot(getDirection())/fSqrLen;
    kDiff -= getDirection()*fT;

    if ( linloc )
      *linloc = fT;

    return kDiff.mag2();
  }
    //stolen from www.geometrictools.com(previously magic-software)
  double Line3::getDistanceSquared(Line3 const &line2,
			      double* pfThisLine, double* pfline2) const
  {

    Vector3d kDiff;
    kDiff = ((*this).getOrigin() - line2.getOrigin());
    double fA01 = -(*this).getDirection().dot(line2.getDirection());
    double fB0 = kDiff.dot((*this).getDirection());
    double fC = kDiff.mag2();
    double fDet = abs((double)1.0 - fA01*fA01);
    double fB1, fS0, fS1, fSqrDist;
    double const ZERO_TOLERANCE = ALMOST_ZERO;
                       //should this be adjustable? fixme
    

    if (fDet >= ZERO_TOLERANCE)
    {
        // lines are not parallel
        fB1 = -kDiff.dot(line2.getDirection());
        double fInvDet = ((double)1.0)/fDet;
        fS0 = (fA01*fB1-fB0)*fInvDet;
        fS1 = (fA01*fB0-fB1)*fInvDet;
        fSqrDist = fS0*(fS0+fA01*fS1+((double)2.0)*fB0) +
            fS1*(fA01*fS0+fS1+((double)2.0)*fB1)+fC;
    }
    else
    {
        // lines are parallel, select any closest pair of points
        fS0 = -fB0;
        fS1 = (double)0.0;
        fSqrDist = fB0*fS0+fC;
    }

    if (pfThisLine)
	*pfThisLine = fS0;
    if (pfline2)
	*pfline2 = fS1;
    return fabs(fSqrDist);
}

  double Line3::getDistanceSquared(Segment3 const &segment,
				   double* pfLinP, double* pfSegP) const
  {
    // alexr, this tolerance is scaled below to avoid
    //        absolute tolerancing issues
    static double DETERMINANT_TOL = 1.0e-9;

    // (www.magic-software.com)
    Vector3d kDiff(getOrigin() - segment.getCoordinate(0));
    double fA00 = getDirection().mag2();
    double fA01 = -getDirection().dot(segment.getDirection());
    double fA11 = segment.getDirection().mag2();
    double fB0 = kDiff.dot(getDirection());
    double fC = kDiff.mag2();
    double fDet = fabs(fA00*fA11-fA01*fA01);
    double fB1, fS, fT, fSqrDist;

    // alexr. 1-16-2013 improving tolerances below, adding first case
    if (fA00 < ALMOST_ZERO)
      {
        // the line is very short so direction is totally unreliable,
        // distance between segment and line origins
        fS = (double)0.0;
        fT = (double)0.0;
        fSqrDist = fC;
      }
    else if (fA11 > ALMOST_ZERO &&  fDet >= DETERMINANT_TOL*fA00*fA11 )
      {
        // case that lines are not parallel: this is the usual intersection check
        fB1 = -kDiff.dot(segment.getDirection());
        fT = fA01*fB0-fA00*fB1;
        
        if ( fT >= (double)0.0 )
          {
            if ( fT <= fDet )
              {
                // two interior points are closest, one on line and one on
                // segment
                double fInvDet = ((double)1.0)/fDet;
                fS = (fA01*fB1-fA11*fB0)*fInvDet;
                fT *= fInvDet;
                fSqrDist = fS*(fA00*fS+fA01*fT+((double)2.0)*fB0) +
                  fT*(fA01*fS+fA11*fT+((double)2.0)*fB1)+fC;
              }
            else
              {
                // end point e1 of segment and interior point of line are
                // closest
                double fTmp = fA01+fB0;
                fS = -fTmp/fA00;
                fT = (double)1.0;
                fSqrDist = fTmp*fS+fA11+((double)2.0)*fB1+fC;
              }
          }
        else
          {
            // end point e0 of segment and interior point of line are closest
            fS = -fB0/fA00;
            fT = (double)0.0;
            fSqrDist = fB0*fS+fC;
          }
      }
    else
      {
        // lines are parallel or the segment is very short,
        // closest pair with one point at segment origin
        fS = -fB0/fA00;
        fT = (double)0.0;
        fSqrDist = fB0*fS+fC;
      }

    if ( pfLinP )
      *pfLinP = fS;
    
    if ( pfSegP )
      *pfSegP = fT;
    
    return fabs(fSqrDist);
  }
}
