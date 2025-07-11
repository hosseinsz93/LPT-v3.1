#include "GeomBox3.h"

#include "GeomConeFrustum3.h"
#include "GeomLine3.h"
#include "GeomRay3.h"
#include "GeomSegment3.h"
#include "GeomSphere3.h"

#include "GeomBoundingBox.h"
#include "GeomUtility.h"

#include <limits>

namespace Geom
{
  Box3::Box3(Vector3d const &corner1,
	     Vector3d const &corner2)
  {
    for (int i=0; i<3; ++i)
      {
	_min(i) = std::min(corner1(i), corner2(i));
	_max(i) = std::max(corner1(i), corner2(i));
      }
  }

  Box3::Box3(std::vector < Vector3d > const &queryFace,
             double adjustment)
  {
    _min(0) = _min(1) = _min(2) =  std::numeric_limits<double>::max();
    _max(0) = _max(1) = _max(2) = -std::numeric_limits<double>::max();

    std::vector < Vector3d >::const_iterator it;
    for (it=queryFace.begin(); it!=queryFace.end(); ++it)
      for (int i=0; i<3; ++i)
      {
        _min(i) = std::min(_min(i), (*it)(i));
        _max(i) = std::max(_max(i), (*it)(i));
      }

    for (int i=0; i<3; ++i)
    {
      _min(i) -= adjustment;
      _max(i) += adjustment;
    }
  }

  Box3::Box3(Box3 const &rhs)
    : BoundedShape3(rhs)
    , _min(rhs._min)
    , _max(rhs._max) 
  {}


  Box3::Box3(BoundingBox const &bbox)
    :_min(bbox.min()),
     _max(bbox.max())
  {}

  BoundingBox Box3::getBoundingBox() const
  {
    return BoundingBox(_min,
                       _max);
  }
  void Box3::inflate(double const size)
  {
    for (int i=0; i<3; i++)
      {
        _min(i) -= size;
        _max(i) += size;
      }
  }

  bool Box3::intersectsBox(BoundingBox const &bbox) const
  {
    return testIntersection(Box3(bbox.min(), bbox.max()));
  }

  bool Box3::containsPoint(Vector3d const &point) const
  {
    return testIntersection(point);
  }

  double Box3::getDistanceSquared(Vector3d const &coord) const
  {
    double distancesquared = 0.0;
    double delta;

    if (coord(0) < this->min()(0))
      {
	delta = coord(0) - this->min()(0);
	distancesquared += delta*delta;
      }
    else if (coord(0) > this->max()(0))
      {
	delta = coord(0) - this->max()(0);
	distancesquared += delta*delta;
      }
    
    if (coord(1) < this->min()(1))
      {
	delta = coord(1) - this->min()(1);
	distancesquared += delta*delta;
      }
    else if (coord(1) > this->max()(1))
      {
	delta = coord(1) - this->max()(1);
	distancesquared += delta*delta;
      }
    
    if (coord(2) < this->min()(2))
      {
	delta = coord(2) - this->min()(2);
	distancesquared += delta*delta;
      }
    else if (coord(2) > this->max()(2))
      {
	delta = coord(2) - this->max()(2);
	distancesquared += delta*delta;
      }
    
    return distancesquared;
  }



  // Stolen almost verbatim from Wild Magic library (www.magic-software.com)
  static bool Clip (double fDenom, double fNumer, double &rfT0, double &rfT1)
  {
    // Return value is 'true' if line segment intersects the current test
    // plane.  Otherwise 'false' is returned in which case the line segment
    // is entirely clipped.
	
    if ( fDenom > 0.0 )
      {
	if ( fNumer > fDenom*rfT1 )
	  return false;
	if ( fNumer > fDenom*rfT0 )
	  rfT0 = fNumer/fDenom;
	return true;
      }
    else if ( fDenom < 0.0 )
      {
	if ( fNumer > fDenom*rfT0 )
	  return false;
	if ( fNumer > fDenom*rfT1 )
	  rfT1 = fNumer/fDenom;
	return true;
      }
    else
      {
	return fNumer <= 0.0;
      }
  }

  // Stolen almost verbatim from Wild Magic library (www.magic-software.com)
  static bool FindIntersection(Vector3d const &origin,
			       Vector3d const &direction,
			       Vector3d const &extent,
			       double &rfT0,
			       double &rfT1)
  {
    double fSaveT0 = rfT0, fSaveT1 = rfT1;
    
    bool bNotEntirelyClipped =
      Clip(+direction(0),-origin(0)-extent(0),rfT0,rfT1) &&
      Clip(-direction(0),+origin(0)-extent(0),rfT0,rfT1) &&
      Clip(+direction(1),-origin(1)-extent(1),rfT0,rfT1) &&
      Clip(-direction(1),+origin(1)-extent(1),rfT0,rfT1) &&
      Clip(+direction(2),-origin(2)-extent(2),rfT0,rfT1) &&
      Clip(-direction(2),+origin(2)-extent(2),rfT0,rfT1);
    
    return bNotEntirelyClipped && ( rfT0 != fSaveT0 || rfT1 != fSaveT1 );
  }

  static inline Box3 createInflatedBox(Box3 const &box,
                                       double inflateFactor)
  {
    Box3 inflateBox(box);
    double dim=inflateBox.getDiagonal().sum()/3.0;
    double inflationAmt = inflateFactor*dim;
    inflateBox.inflate(inflationAmt);

    return inflateBox;
  }
  bool Box3::testIntersection(Segment3 const &segment,
                              double inflateFactor) const
  {	
    //NOTE:  we will use inflateFactor to inflate the bbox by a 
    //small amount to ensure that
    //segments passing through a corner edge of the box will be
    //found as intersecting.  We saw problems with octree refinement
    //at segments because especially as the box shrank, we saw roundoff
    //which caused segment to no longer be intersecting
    Box3 inflateBox;
    Box3 const *testBox;

    if (inflateFactor != 0.0)
    {
      inflateBox = createInflatedBox(*this,inflateFactor);
      testBox = &inflateBox;
    }
    else
      testBox = this;

    //------------------------------------------------------------------------------------------
    // If bbox of segment is not intersecting return false.
    
    Box3 tempBox(segment.getConstBoundingBox());
    if (!testBox->testIntersection(tempBox))
      {
	return false;
      }
   
 
    if (testBox->testIntersection(segment.getCoordinate(0)) ||
	testBox->testIntersection(segment.getCoordinate(1)))
    {
      return true;
    }
    
    Vector3d extent = testBox->getDiagonal();
    extent *= 0.5;

    Vector3d kDiff(segment.getCoordinate(0) - testBox->getCenter());

    double fT0 = 0.0;
    double fT1 = segment.getLength();
    if (!FindIntersection(kDiff, segment.getUnitDirection(), extent, 
			  fT0, fT1))
      return false;
    
    return true;
  }
			
  static void Face (int i0, int i1, int i2, Vector3d &rkPnt,
		    Vector3d const &rkDir, Box3 const &rkBox,
		    Vector3d const &rkPmE, double* pfLParam, 
		    double &rfSqrDistance)
  {
    Vector3d kPpE;
    double fLSqr, fInv, fTmp, fParam, fT, fDelta;

    

    kPpE(i1) = rkPnt(i1) + rkBox.getExtent(i1);
    kPpE(i2) = rkPnt(i2) + rkBox.getExtent(i2);
    if ( rkDir(i0)*kPpE(i1) >= rkDir(i1)*rkPmE(i0) )
      {
        if ( rkDir(i0)*kPpE(i2) >= rkDir(i2)*rkPmE(i0) )
	  {
            // v[i1] >= -e[i1], v[i2] >= -e[i2] (distance = 0)
            if ( pfLParam )
	      {
                rkPnt(i0) = rkBox.getExtent(i0);
                fInv = 1.0/rkDir(i0);
                rkPnt(i1) -= rkDir(i1)*rkPmE(i0)*fInv;
                rkPnt(i2) -= rkDir(i2)*rkPmE(i0)*fInv;
                *pfLParam = -rkPmE(i0)*fInv;
	      }
	  }
        else
	  {
            // v[i1] >= -e[i1], v[i2] < -e[i2]
            fLSqr = rkDir(i0)*rkDir(i0) + rkDir(i2)*rkDir(i2);
            fTmp = fLSqr*kPpE(i1) - rkDir(i1)*(rkDir(i0)*rkPmE(i0) +
					       rkDir(i2)*kPpE(i2));
            if ( fTmp <= 2.0*fLSqr*rkBox.getExtent(i1) )
	      {
                fT = fTmp/fLSqr;
                fLSqr += rkDir(i1)*rkDir(i1);
                fTmp = kPpE(i1) - fT;
                fDelta = rkDir(i0)*rkPmE(i0) + rkDir(i1)*fTmp +
		  rkDir(i2)*kPpE(i2);
                fParam = -fDelta/fLSqr;
                rfSqrDistance += rkPmE(i0)*rkPmE(i0) + fTmp*fTmp +
		  kPpE(i2)*kPpE(i2) + fDelta*fParam;
		
                if ( pfLParam )
		  {
                    *pfLParam = fParam;
                    rkPnt(i0) = rkBox.getExtent(i0);
                    rkPnt(i1) = fT - rkBox.getExtent(i1);
                    rkPnt(i2) = -rkBox.getExtent(i2);
		  }
	      }
            else
	      {
                fLSqr += rkDir(i1)*rkDir(i1);
                fDelta = rkDir(i0)*rkPmE(i0) + rkDir(i1)*rkPmE(i1) +
		  rkDir(i2)*kPpE(i2);
                fParam = -fDelta/fLSqr;
                rfSqrDistance += rkPmE(i0)*rkPmE(i0) + rkPmE(i1)*rkPmE(i1) +
		  kPpE(i2)*kPpE(i2) + fDelta*fParam;
		
                if ( pfLParam )
		  {
                    *pfLParam = fParam;
                    rkPnt(i0) = rkBox.getExtent(i0);
                    rkPnt(i1) = rkBox.getExtent(i1);
                    rkPnt(i2) = -rkBox.getExtent(i2);
		  }
	      }
	  }
      }
    else
      {
        if ( rkDir(i0)*kPpE(i2) >= rkDir(i2)*rkPmE(i0) )
	  {
            // v[i1] < -e[i1], v[i2] >= -e[i2]
            fLSqr = rkDir(i0)*rkDir(i0) + rkDir(i1)*rkDir(i1);
            fTmp = fLSqr*kPpE(i2) - rkDir(i2)*(rkDir(i0)*rkPmE(i0) +
					       rkDir(i1)*kPpE(i1));
            if ( fTmp <= 2.0*fLSqr*rkBox.getExtent(i2) )
	      {
                fT = fTmp/fLSqr;
                fLSqr += rkDir(i2)*rkDir(i2);
                fTmp = kPpE(i2) - fT;
                fDelta = rkDir(i0)*rkPmE(i0) + rkDir(i1)*kPpE(i1) +
		  rkDir(i2)*fTmp;
                fParam = -fDelta/fLSqr;
                rfSqrDistance += rkPmE(i0)*rkPmE(i0) + kPpE(i1)*kPpE(i1) +
		  fTmp*fTmp + fDelta*fParam;
		
                if ( pfLParam )
		  {
                    *pfLParam = fParam;
                    rkPnt(i0) = rkBox.getExtent(i0);
                    rkPnt(i1) = -rkBox.getExtent(i1);
                    rkPnt(i2) = fT - rkBox.getExtent(i2);
		  }
	      }
            else
	      {
                fLSqr += rkDir(i2)*rkDir(i2);
                fDelta = rkDir(i0)*rkPmE(i0) + rkDir(i1)*kPpE(i1) +
		  rkDir(i2)*rkPmE(i2);
                fParam = -fDelta/fLSqr;
                rfSqrDistance += rkPmE(i0)*rkPmE(i0) + kPpE(i1)*kPpE(i1) +
		  rkPmE(i2)*rkPmE(i2) + fDelta*fParam;
		
                if ( pfLParam )
		  {
                    *pfLParam = fParam;
                    rkPnt(i0) = rkBox.getExtent(i0);
                    rkPnt(i1) = -rkBox.getExtent(i1);
                    rkPnt(i2) = rkBox.getExtent(i2);
		  }
	      }
	  }
        else
	  {
            // v[i1] < -e[i1], v[i2] < -e[i2]
            fLSqr = rkDir(i0)*rkDir(i0)+rkDir(i2)*rkDir(i2);
            fTmp = fLSqr*kPpE(i1) - rkDir(i1)*(rkDir(i0)*rkPmE(i0) +
					       rkDir(i2)*kPpE(i2));
            if ( fTmp >= 0.0 )
	      {
                // v[i1]-edge is closest
                if ( fTmp <= 2.0*fLSqr*rkBox.getExtent(i1) )
		  {
                    fT = fTmp/fLSqr;
                    fLSqr += rkDir(i1)*rkDir(i1);
                    fTmp = kPpE(i1) - fT;
                    fDelta = rkDir(i0)*rkPmE(i0) + rkDir(i1)*fTmp +
		      rkDir(i2)*kPpE(i2);
                    fParam = -fDelta/fLSqr;
                    rfSqrDistance += rkPmE(i0)*rkPmE(i0) + fTmp*fTmp +
		      kPpE(i2)*kPpE(i2) + fDelta*fParam;
		    
                    if ( pfLParam )
		      {
                        *pfLParam = fParam;
                        rkPnt(i0) = rkBox.getExtent(i0);
                        rkPnt(i1) = fT - rkBox.getExtent(i1);
                        rkPnt(i2) = -rkBox.getExtent(i2);
		      }
		  }
                else
		  {
                    fLSqr += rkDir(i1)*rkDir(i1);
                    fDelta = rkDir(i0)*rkPmE(i0) + rkDir(i1)*rkPmE(i1) +
		      rkDir(i2)*kPpE(i2);
                    fParam = -fDelta/fLSqr;
                    rfSqrDistance += rkPmE(i0)*rkPmE(i0) + rkPmE(i1)*rkPmE(i1)
		      + kPpE(i2)*kPpE(i2) + fDelta*fParam;
		    
                    if ( pfLParam )
		      {
                        *pfLParam = fParam;
                        rkPnt(i0) = rkBox.getExtent(i0);
                        rkPnt(i1) = rkBox.getExtent(i1);
                        rkPnt(i2) = -rkBox.getExtent(i2);
		      }
		  }
                return;
	      }
	    
            fLSqr = rkDir(i0)*rkDir(i0) + rkDir(i1)*rkDir(i1);
            fTmp = fLSqr*kPpE(i2) - rkDir(i2)*(rkDir(i0)*rkPmE(i0) +
					       rkDir(i1)*kPpE(i1));
            if ( fTmp >= 0.0 )
	      {
                // v[i2]-edge is closest
                if ( fTmp <= 2.0*fLSqr*rkBox.getExtent(i2) )
		  {
                    fT = fTmp/fLSqr;
                    fLSqr += rkDir(i2)*rkDir(i2);
                    fTmp = kPpE(i2) - fT;
                    fDelta = rkDir(i0)*rkPmE(i0) + rkDir(i1)*kPpE(i1) +
		      rkDir(i2)*fTmp;
                    fParam = -fDelta/fLSqr;
                    rfSqrDistance += rkPmE(i0)*rkPmE(i0) + kPpE(i1)*kPpE(i1) +
		      fTmp*fTmp + fDelta*fParam;
		    
                    if ( pfLParam )
		      {
                        *pfLParam = fParam;
                        rkPnt(i0) = rkBox.getExtent(i0);
                        rkPnt(i1) = -rkBox.getExtent(i1);
                        rkPnt(i2) = fT - rkBox.getExtent(i2);
		      }
		  }
                else
		  {
                    fLSqr += rkDir(i2)*rkDir(i2);
                    fDelta = rkDir(i0)*rkPmE(i0) + rkDir(i1)*kPpE(i1) +
		      rkDir(i2)*rkPmE(i2);
                    fParam = -fDelta/fLSqr;
                    rfSqrDistance += rkPmE(i0)*rkPmE(i0) + kPpE(i1)*kPpE(i1) +
		      rkPmE(i2)*rkPmE(i2) + fDelta*fParam;
		    
                    if ( pfLParam )
		      {
                        *pfLParam = fParam;
                        rkPnt(i0) = rkBox.getExtent(i0);
                        rkPnt(i1) = -rkBox.getExtent(i1);
                        rkPnt(i2) = rkBox.getExtent(i2);
		      }
		  }
                return;
	      }
	    
            // (v[i1],v[i2])-corner is closest
            fLSqr += rkDir(i2)*rkDir(i2);
            fDelta = rkDir(i0)*rkPmE(i0) + rkDir(i1)*kPpE(i1) +
	      rkDir(i2)*kPpE(i2);
            fParam = -fDelta/fLSqr;
            rfSqrDistance += rkPmE(i0)*rkPmE(i0) + kPpE(i1)*kPpE(i1) +
	      kPpE(i2)*kPpE(i2) + fDelta*fParam;
	    
            if ( pfLParam )
	      {
                *pfLParam = fParam;
                rkPnt(i0) = rkBox.getExtent(i0);
                rkPnt(i1) = -rkBox.getExtent(i1);
                rkPnt(i2) = -rkBox.getExtent(i2);
	      }
	  }
      }
  }

  static void CaseNoZeros (Vector3d &rkPnt, Vector3d const &rkDir,
			   Box3 const &rkBox, double *pfLParam, 
			   double &rfSqrDistance)
  {
    Vector3d extent = rkBox.getDiagonal();
    extent *= 0.5;

    Vector3d kPmE(rkPnt - extent);

    double fProdDxPy = rkDir(0)*kPmE(1);
    double fProdDyPx = rkDir(1)*kPmE(0);
    double fProdDzPx, fProdDxPz, fProdDzPy, fProdDyPz;
    
    if ( fProdDyPx >= fProdDxPy )
      {
        fProdDzPx = rkDir(2)*kPmE(0);
        fProdDxPz = rkDir(0)*kPmE(2);
        if ( fProdDzPx >= fProdDxPz )
	  {
            // line intersects x = e0
            Face(0,1,2,rkPnt,rkDir,rkBox,kPmE,pfLParam,rfSqrDistance);
	  }
        else
	  {
            // line intersects z = e2
            Face(2,0,1,rkPnt,rkDir,rkBox,kPmE,pfLParam,rfSqrDistance);
	  }
      }
    else
      {
        fProdDzPy = rkDir(2)*kPmE(1);
        fProdDyPz = rkDir(1)*kPmE(2);
        if ( fProdDzPy >= fProdDyPz )
	  {
            // line intersects y = e1
            Face(1,2,0,rkPnt,rkDir,rkBox,kPmE,pfLParam,rfSqrDistance);
	  }
        else
	  {
            // line intersects z = e2
            Face(2,0,1,rkPnt,rkDir,rkBox,kPmE,pfLParam,rfSqrDistance);
	  }
      }
  }


  static void Case0 (int i0, int i1, int i2, Vector3d &rkPnt,
		     Vector3d const &rkDir, Box3 const &rkBox, 
		     double *pfLParam, double &rfSqrDistance)
  {
    double fPmE0 = rkPnt(i0) - rkBox.getExtent(i0);
    double fPmE1 = rkPnt(i1) - rkBox.getExtent(i1);
    double fProd0 = rkDir(i1)*fPmE0;
    double fProd1 = rkDir(i0)*fPmE1;
    double fDelta, fInvLSqr, fInv;

    if ( fProd0 >= fProd1 )
      {
        // line intersects P[i0] = e[i0]
        rkPnt(i0) = rkBox.getExtent(i0);
	
        double fPpE1 = rkPnt(i1) + rkBox.getExtent(i1);
        fDelta = fProd0 - rkDir(i0)*fPpE1;
        if ( fDelta >= 0.0 )
	  {
            fInvLSqr = 1.0/(rkDir(i0)*rkDir(i0)+rkDir(i1)*rkDir(i1));
            rfSqrDistance += fDelta*fDelta*fInvLSqr;
            if ( pfLParam )
	      {
                rkPnt(i1) = -rkBox.getExtent(i1);
                *pfLParam = -(rkDir(i0)*fPmE0+rkDir(i1)*fPpE1)*fInvLSqr;
	      }
	  }
        else
	  {
            if ( pfLParam )
	      {
                fInv = 1.0/rkDir(i0);
                rkPnt(i1) -= fProd0*fInv;
                *pfLParam = -fPmE0*fInv;
	      }
	  }
      }
    else
      {
        // line intersects P[i1] = e[i1]
        rkPnt(i1) = rkBox.getExtent(i1);

        double fPpE0 = rkPnt(i0) + rkBox.getExtent(i0);
        fDelta = fProd1 - rkDir(i1)*fPpE0;
        if ( fDelta >= 0.0 )
	  {
            fInvLSqr = 1.0/(rkDir(i0)*rkDir(i0)+rkDir(i1)*rkDir(i1));
            rfSqrDistance += fDelta*fDelta*fInvLSqr;
            if ( pfLParam )
	      {
                rkPnt(i0) = -rkBox.getExtent(i0);
                *pfLParam = -(rkDir(i0)*fPpE0+rkDir(i1)*fPmE1)*fInvLSqr;
	      }
	  }
        else
	  {
            if ( pfLParam )
	      {
                fInv = 1.0/rkDir(i1);
                rkPnt(i0) -= fProd1*fInv;
                *pfLParam = -fPmE1*fInv;
	      }
	  }
      }
    
    if ( rkPnt(i2) < -rkBox.getExtent(i2))
      {
        fDelta = rkPnt(i2) + rkBox.getExtent(i2);
        rfSqrDistance += fDelta*fDelta;
        rkPnt(i2) = -rkBox.getExtent(i2);
      }
    else if ( rkPnt(i2) > rkBox.getExtent(i2) )
      {
        fDelta = rkPnt(i2) - rkBox.getExtent(i2);
        rfSqrDistance += fDelta*fDelta;
        rkPnt(i2) = rkBox.getExtent(i2);
      }
  }
  
  static void Case00 (int i0, int i1, int i2, Vector3d &rkPnt,
		      Vector3d const &rkDir, Box3 const &rkBox, 
		      double *pfLParam, double &rfSqrDistance)
  {
    double fDelta;
    
    if ( pfLParam )
      *pfLParam = (rkBox.getExtent(i0) - rkPnt(i0))/rkDir(i0);
    
    rkPnt(i0) = rkBox.getExtent(i0);
    
    if ( rkPnt(i1) < -rkBox.getExtent(i1) )
      {
        fDelta = rkPnt(i1) + rkBox.getExtent(i1);
	rfSqrDistance += fDelta*fDelta;
	rkPnt(i1) = -rkBox.getExtent(i1);
      }
    else if ( rkPnt(i1) > rkBox.getExtent(i1) )
      {
	fDelta = rkPnt(i1) - rkBox.getExtent(i1);
        rfSqrDistance += fDelta*fDelta;
        rkPnt(i1) = rkBox.getExtent(i1);
      }
    
    if ( rkPnt(i2) < -rkBox.getExtent(i2) )
      {
        fDelta = rkPnt(i2) + rkBox.getExtent(i2);
        rfSqrDistance += fDelta*fDelta;
        rkPnt(i1) = -rkBox.getExtent(i2);
      }
    else if ( rkPnt(i2) > rkBox.getExtent(i2) )
      {
        fDelta = rkPnt(i2) - rkBox.getExtent(i2);
        rfSqrDistance += fDelta*fDelta;
        rkPnt(i2) = rkBox.getExtent(i2);
      }
  }

  static void Case000 (Vector3d &rkPnt, Box3 const &rkBox,
		       double &rfSqrDistance)
  {
    double fDelta;
    
    if ( rkPnt(0) < -rkBox.getExtent(0) )
      {
        fDelta = rkPnt(0) + rkBox.getExtent(0);
        rfSqrDistance += fDelta*fDelta;
        rkPnt(0) = -rkBox.getExtent(0);
      }
    else if ( rkPnt(0) > rkBox.getExtent(0) )
      {
        fDelta = rkPnt(0) - rkBox.getExtent(0);
        rfSqrDistance += fDelta*fDelta;
        rkPnt(0) = rkBox.getExtent(0);
      }
    
    if ( rkPnt(1) < -rkBox.getExtent(1) )
      {
        fDelta = rkPnt(1) + rkBox.getExtent(1);
        rfSqrDistance += fDelta*fDelta;
        rkPnt(1) = -rkBox.getExtent(1);
      }
    else if ( rkPnt(1) > rkBox.getExtent(1) )
      {
        fDelta = rkPnt(1) - rkBox.getExtent(1);
        rfSqrDistance += fDelta*fDelta;
        rkPnt(1) = rkBox.getExtent(1);
      }
    
    if ( rkPnt(2) < -rkBox.getExtent(2) )
      {
        fDelta = rkPnt(2) + rkBox.getExtent(2);
        rfSqrDistance += fDelta*fDelta;
        rkPnt(2) = -rkBox.getExtent(2);
      }
    else if ( rkPnt(2) > rkBox.getExtent(2) )
      {
        fDelta = rkPnt(2) - rkBox.getExtent(2);
        rfSqrDistance += fDelta*fDelta;
        rkPnt(2) = rkBox.getExtent(2);
      }
  }

  static double SqrDistance (Line3 const &rkLine, Box3 const &rkBox, double *pfLParam,
			     double *pfBParam0, double *pfBParam1, double *pfBParam2)
  {
    // compute coordinates of line in box coordinate system
    Vector3d kDiff(rkLine.getOrigin() - rkBox.getCenter());
    Vector3d kPnt = kDiff;
    Vector3d kDir = rkLine.getDirection();

    // Apply reflections so that direction vector has nonnegative components.
    bool bReflect[3];
    for (unsigned i = 0; i < 3; i++)
      {
        if ( kDir(i) < 0.0 )
	  {
            kPnt(i) = -kPnt(i);
            kDir(i) = -kDir(i);
            bReflect[i] = true;
	  }
        else
	  {
            bReflect[i] = false;
	  }
      }

    double fSqrDistance = 0.0;
    
    if ( kDir(0) > 0.0 )
      { 
        if ( kDir(1) > 0.0 )
	  {
            if ( kDir(2) > 0.0 )
	      {
                // (+,+,+)
                CaseNoZeros(kPnt,kDir,rkBox,pfLParam,fSqrDistance);
	      }
            else
	      {
                // (+,+,0)
                Case0(0,1,2,kPnt,kDir,rkBox,pfLParam,fSqrDistance);
	      }
	  }
        else
	  {
            if ( kDir(2) > 0.0 )
	      {
                // (+,0,+)
                Case0(0,2,1,kPnt,kDir,rkBox,pfLParam,fSqrDistance);
	      }
            else
	      {
                // (+,0,0)
                Case00(0,1,2,kPnt,kDir,rkBox,pfLParam,fSqrDistance);
	      }
	  }
      }
    else
      {
        if ( kDir(1) > 0.0 )
	  {
            if ( kDir(2) > 0.0 )
	      {
                // (0,+,+)
                Case0(1,2,0,kPnt,kDir,rkBox,pfLParam,fSqrDistance);
	      }
            else
	      {
                // (0,+,0)
                Case00(1,0,2,kPnt,kDir,rkBox,pfLParam,fSqrDistance);
	      }
	  }
        else
	  {
            if ( kDir(2) > 0.0 )
	      {
                // (0,0,+)
                Case00(2,0,1,kPnt,kDir,rkBox,pfLParam,fSqrDistance);
	      }
            else
	      {
                // (0,0,0)
                Case000(kPnt,rkBox,fSqrDistance);
                if ( pfLParam )
		  *pfLParam = 0.0;
	      }
	  }
      }

    // undo reflections
    for (unsigned i = 0; i < 3; i++)
      {
        if ( bReflect[i] )
	  kPnt(i) = -kPnt(i);
      }
    
    if ( pfBParam0 )
      *pfBParam0 = kPnt(0);

    if ( pfBParam1 )
      *pfBParam1 = kPnt(1);
    
    if ( pfBParam2 )
      *pfBParam2 = kPnt(2);
    
    return fSqrDistance;
  }


  static double SqrDistance(Vector3d const &rkPoint, Box3 const &rkBox,
			    double *pfBParam0, double *pfBParam1, double *pfBParam2)
  {
    // compute coordinates of point in box coordinate system
    Vector3d kDiff(rkPoint - rkBox.getCenter());
    Vector3d kClosest = kDiff;
    
    // project test point onto box
    double fSqrDistance = 0.0;
    double fDelta;
 
    if ( kClosest(0) < -rkBox.getExtent(0) )
      {
        fDelta = kClosest(0) + rkBox.getExtent(0);
        fSqrDistance += fDelta*fDelta;
        kClosest(0) = -rkBox.getExtent(0);
      }
    else if ( kClosest(0) > rkBox.getExtent(0) )
      {
        fDelta = kClosest(0) - rkBox.getExtent(0);
        fSqrDistance += fDelta*fDelta;
        kClosest(0) = rkBox.getExtent(0);
      }
    
    if ( kClosest(1) < -rkBox.getExtent(1) )
      {
        fDelta = kClosest(1) + rkBox.getExtent(1);
        fSqrDistance += fDelta*fDelta;
        kClosest(1) = -rkBox.getExtent(1);
      }
    else if ( kClosest(1) > rkBox.getExtent(1) )
      {
        fDelta = kClosest(1) - rkBox.getExtent(1);
        fSqrDistance += fDelta*fDelta;
        kClosest(1) = rkBox.getExtent(1);
      }
    
    if ( kClosest(2) < -rkBox.getExtent(2) )
      {
        fDelta = kClosest(2) + rkBox.getExtent(2);
        fSqrDistance += fDelta*fDelta;
        kClosest(2) = -rkBox.getExtent(2);
      }
    else if ( kClosest(2) > rkBox.getExtent(2) )
      {
        fDelta = kClosest(2) - rkBox.getExtent(2);
        fSqrDistance += fDelta*fDelta;
        kClosest(2) = rkBox.getExtent(2);
      }
    
    if ( pfBParam0 )
      *pfBParam0 = kClosest(0);
    
    if ( pfBParam1 )
      *pfBParam1 = kClosest(1);
                             
    if ( pfBParam2 )
      *pfBParam2 = kClosest(2);
    
    return fSqrDistance;
  }

  static double SqrDistance (Segment3 const &rkSeg, Box3 const &rkBox, 
			     double *pfLParam, double *pfBParam0, double *pfBParam1, 
			     double *pfBParam2)
  {
    Line3 kLine(rkSeg.getCoordinate(0),
		Vector3d(rkSeg.getCoordinate(1) - rkSeg.getCoordinate(0)));
    
    double fLP, fBP0, fBP1, fBP2;
    double fSqrDistance = SqrDistance(kLine,rkBox,&fLP,&fBP0,&fBP1,&fBP2);
    if ( fLP >= 0.0 )
      {
        if ( fLP <= 1.0 )
	  {
            if ( pfLParam )
	      *pfLParam = fLP;
	    
            if ( pfBParam0 )
	      *pfBParam0 = fBP0;
	    
            if ( pfBParam1 )
	      *pfBParam1 = fBP1;
	    
            if ( pfBParam2 )
	      *pfBParam2 = fBP2;
	    
            return fSqrDistance;
	  }
        else
	  {
	    Vector3d vec = rkSeg.getCoordinate(1);
            fSqrDistance = SqrDistance(vec,rkBox,pfBParam0,pfBParam1,
				       pfBParam2);
	    
            if ( pfLParam )
	      *pfLParam = 1.0;
	    
            return fSqrDistance;
	  }
      }
    else
      {
	Vector3d vec = rkSeg.getCoordinate(0);
        fSqrDistance = SqrDistance(vec,rkBox,pfBParam0,pfBParam1,pfBParam2);
	
        if ( pfLParam )
	  *pfLParam = 0.0;
	
        return fSqrDistance;
      }
  }

  double Box3::getDistanceSquared(Segment3 const &segment3) const
  {
    double fLP, fBP0, fBP1, fBP2;
    double d = SqrDistance(segment3,*this,&fLP,&fBP0,&fBP1,&fBP2);
    // the above algorithm can sometimes return very small, negative numbers
    // when it should return 0.0.  Since the return value is often passed to
    // sqrt, make sure it's not negative so we don't get a FPE.  Note, this
    // needs to be revisited later to make sure we're not getting any
    // wrong answers here!
    return std::max(0.0, d);
  }


  bool Box3::testIntersection(Sphere3 const &sphere) const
  {
    Vector3d boxcenter = this->getCenter();
    Vector3d boxextent = this->getDiagonal();
    boxextent *= 0.5;

    double fAx = abs(sphere.getCenter()(0) - boxcenter(0));
    double fAy = abs(sphere.getCenter()(1) - boxcenter(1));
    double fAz = abs(sphere.getCenter()(2) - boxcenter(2));
    double fDx = fAx - boxextent(0);
    double fDy = fAy - boxextent(1);
    double fDz = fAz - boxextent(2);

    double radius = sphere.getRadius();

    if ( fAx <= boxextent(0) )
      {
        if ( fAy <= boxextent(1) )
          {
            if ( fAz <= boxextent(2) )
              {
                // sphere center inside box
                return true;
              }
            else
              {
                // potential sphere-face intersection with face z
                return fDz <= radius;
              }
          }
        else
          {
            if ( fAz <= boxextent(2) )
              {
                // potential sphere-face intersection with face y
                return fDy <= radius;
              }
            else
              {
                // potential sphere-edge intersection with edge formed
                // by faces y and z
                double fRSqr = radius*radius;
                return fDy*fDy + fDz*fDz <= fRSqr;
              }
          }
      }
    else
      {
        if ( fAy <= boxextent(1) )
          {
            if ( fAz <= boxextent(2) )
              {
                // potential sphere-face intersection with face x
                return fDx <= radius;
              }
            else
              {
                // potential sphere-edge intersection with edge formed
                // by faces x and z
                double fRSqr = radius*radius;
                return fDx*fDx + fDz*fDz <= fRSqr;
              }
          }
        else
          {
            if ( fAz <= boxextent(2) )
              {
                // potential sphere-edge intersection with edge formed
                // by faces x and y
                double fRSqr = radius*radius;
                return fDx*fDx + fDy*fDy <= fRSqr;
              }
            else
              {
                // potential sphere-vertex intersection at corner formed
                // by faces x,y,z
                double fRSqr = radius*radius;
                return fDx*fDx + fDy*fDy + fDz*fDz <= fRSqr;
              }
          }
      }       
    // can't get here...
    return false;
  }

  /*
     * WARNING!!!
     * This test is not completely precise.  It may get some bboxes that
     * are right around one end of the frustum but just outside the frustum
     * end cap.  This issue has been in here from the beginning so I am
     * not going to fix it the day before the cut 4.02 as I don't have a fix
     * worked out and it needs sufficient testing.  For the purposes of
     * refinement, it is fine but not precise and needs to be fixed.  
     *  I have assigned myself a bug report CCMP-6333
     * Chris Geisert 01/05/2008
     */

  bool Box3::testIntersection(ConeFrustum3 const &frustum) const
  {
    enum {
      kPtOutsideConeWithinAxis,
      kPtInsideConeBelowAxis,
      kPtInsideConeAboveAxis,
      kPtOutsideConeBelowAxis,
      kPtOutsideConeAboveAxis,
      kPtBelowConeVertex
    };
    typedef struct _ptInfo {
        Vector3d pt;
      int status;
      bool onBox;
    } ptInfo;
    ptInfo ptList[8];
    
    BoundingBox frustumBox = frustum.getBoundingBox();
    
    if (!Box3(frustumBox.min(), frustumBox.max()).testIntersection(*this))
      return false;
    //find intersection of the two boxes.
    Vector3d boxMax,boxMin,frustBox[2],newBox[2];
    
    boxMax = this->max();
    boxMin = this->min();
    frustBox[1] = frustumBox.max();
    frustBox[0] = frustumBox.min();
    
    newBox[0] = Vector3d(std::max(frustBox[0](0),boxMin(0)),
				 std::max(frustBox[0](1),boxMin(1)),
				 std::max(frustBox[0](2),boxMin(2)));
    newBox[1] = Vector3d(std::min(frustBox[1](0),boxMax(0)),
				 std::min(frustBox[1](1),boxMax(1)),
				 std::min(frustBox[1](2),boxMax(2)));
    
    //frustum's box is entirely contained within bbox so must
    //intersect
    if (Vector3d(newBox[1] - frustBox[1]).norm1() < ALMOST_ZERO &&
	Vector3d(newBox[0] - frustBox[0]).norm1() < ALMOST_ZERO)
      return true;
    
    //first we see if any point on the bbox intersects frustum and
    //classify the points if they don't intersect
    int numInsideConeAboveAxis,numInsideConeBelowAxis,
      numOutsideConeAboveAxis,numOutsideConeBelowAxis,
      numOutsideConeWithinAxis,numBelowVertex=0,numOnBBox;
    int ptIndex=0;
    numOnBBox = 0;
    numInsideConeAboveAxis = 0;
    numInsideConeBelowAxis = 0;
    numOutsideConeAboveAxis = 0;
    numOutsideConeBelowAxis = 0;
    numOutsideConeWithinAxis = 0;
    
    Vector3d frustumAxis = frustum.getAxisVector();

    double frustumAxisLength = frustumAxis.normalize();

    double r1 = frustum.getRadius1();
    double r2 = frustum.getRadius2();
    Vector3d const &center1 = frustum.getCenter1();
    Vector3d const &center2 = frustum.getCenter2();
    double tanTheta = (r2-r1)/frustumAxisLength;
    for (int i=0;i<2;i++)
      {
	for (int j=0;j<2;j++)
	  {
	    for (int k=0;k<2;k++,ptIndex++)
	      {
		double vecMagSquared,dotproduct;
		Vector3d pt(newBox[i](0),newBox[j](1),newBox[k](2));
		ptList[ptIndex].pt = pt;
		Vector3d vec(pt - center1);
		vecMagSquared = vec.mag2();
		dotproduct = vec.dot(frustumAxis);
		double radius = dotproduct*tanTheta + r1;
		double radiusOfPointSq = 
		  vecMagSquared - dotproduct*dotproduct;
		if (dotproduct >= 0 && 
		   dotproduct <= frustumAxisLength)//inside caps of frustum
		  {
		    if ( radiusOfPointSq <= radius*radius)
		      {//inside radius of frustum
			return true;
		      }
		    ptList[ptIndex].status = kPtOutsideConeWithinAxis;
		    ++numOutsideConeWithinAxis;
		  }
		else if ((radius >= 0) && (radiusOfPointSq <= radius*radius))
		  {
		    if (dotproduct < 0)//below frustum cap
		      {
			ptList[ptIndex].status = kPtInsideConeBelowAxis;
			++numInsideConeBelowAxis;
		      }
		    else if (dotproduct > frustumAxisLength)//above frustum
		      {
			ptList[ptIndex].status = kPtInsideConeAboveAxis;
			++numInsideConeAboveAxis;
		      }
		  }
		else 
		  {
		    if (dotproduct < 0)
		      {
			if (radius < 0.0)
			  {
			    ptList[ptIndex].status = kPtBelowConeVertex;
			    ++numBelowVertex;
			  }
			else
			  {
			    ptList[ptIndex].status = 
			      kPtOutsideConeBelowAxis;
			    ++numOutsideConeBelowAxis;
			  }
		      }
		    else if (dotproduct > frustumAxisLength)
		      {
			ptList[ptIndex].status = kPtOutsideConeAboveAxis;
			++numOutsideConeAboveAxis;
		      }
		  }
		if (pt(0) == frustBox[0](0) ||
		   pt(1) == frustBox[0](1) ||
		   pt(2) == frustBox[0](2) ||
		   pt(0) == frustBox[1](0) ||
		   pt(1) == frustBox[1](1) ||
		   pt(2) == frustBox[1](2))
		  {
		    ptList[ptIndex].onBox = true;
		    ++numOnBBox;
		  }
		else
		  ptList[ptIndex].onBox = false;
 	      }
	  }
      }
    //all above frustum cap or below frustum cap  then can't intersect
    if (((numInsideConeAboveAxis + numOutsideConeAboveAxis) == 8) || 
       ((numInsideConeBelowAxis + numOutsideConeBelowAxis + 
	 numBelowVertex)  == 8))
      return false;

    //if axis of cone intersects box then must intersect.
    Segment3 axisSeg(frustum.getCenter1(),frustum.getCenter2());
    if (this->testIntersection(axisSeg))
      return true;
    /*		    if(numOutsideConeWithinAxis  && 
      (numInsideConeAboveAxis || numInsideConeBelowAxis))
      return true;
    */
    //if we get here and we have any points above or below frustum caps but
    //inside the infinite cone then some part of bbox must intersect
    if (numInsideConeAboveAxis || numInsideConeBelowAxis)
      return true;
    
    int v0,v1;
    double dist;
    int const edgeArray[12][2] = {{0,1},{2,3},{4,5},{6,7},
				  {0,2},{1,3},{4,6},{5,7},
				  {0,4},{1,5},{2,6},{3,7}};

    //intersect edges - treat r2==r1 as a cylinder

    //Frustum is a CYLINDER
    if (abs(r2-r1) < ALMOST_ZERO)
      { 
	//now look at every edge
	//at this point if any edge intersects the infinite cylinder then we have to
	//intersect - draw cone bbox and you will see
	
	Segment3 axisSeg(frustum.getCenter1(),frustum.getCenter2());
	Vector3d c1, c2, pt, vec;
	double dotproduct,ptRadiusSq,t1,t2;
	double radiusSq = r1*r1;
	for (int e=0;e<12;e++)
	  {
	    v0 = edgeArray[e][0];
	    v1 = edgeArray[e][1];
	    Vector3d vec0 = ptList[v0].pt;
	    Vector3d vec1 = ptList[v1].pt;
	    Segment3 edgeSeg(vec0,vec1);
	    
	    edgeSeg.getDistanceSquared(axisSeg,&t1,&t2);
	    c1 = edgeSeg.getCoordinate(0);
	    c2 = edgeSeg.getCoordinate(1);
	    pt = c1 + (c2-c1)*t1;
	    
	    vec = pt - center1;
	    dotproduct = vec.dot(frustumAxis);
	    
	    ptRadiusSq = vec.mag2() - dotproduct*dotproduct;
	    if (ptRadiusSq < radiusSq)
	      {
		return true;
	      }
	    
	  }
	//at this point if the bbox intersects any part of the frustum axis
	//then it must intersect this is the final weird case to check
	// 
	BoundingBox newBBox(newBox[0],newBox[1]);
	if (newBBox.containsPoint(center1) || 
	   newBBox.containsPoint(center2) ||
	   newBBox.intersectsRay(center1,frustumAxis,&dist)||
	   //need intersectsLine method
            newBBox.intersectsRay(center1,Vector3d(-frustumAxis),&dist))
	  return true;
      }
    else//not a cylinder
      {
	//now look at every edge
	//at this point if any edge intersects the infinite cone then we have to
	//intersect - draw cone bbox and you will see
	double AdotE,magESquared,c0,c1,c2,EdotDelta,deltaDotA;
	double vertexDist,cosThetaSqr;
	int status0,status1;
	Vector3d E, delta, coneVertex;
	
	vertexDist = r1/tanTheta;
	coneVertex = center1 - frustumAxis*vertexDist;
	cosThetaSqr = 1/(1+tanTheta*tanTheta);
	for (int e=0;e<12;e++)
	  {
	    v0 = edgeArray[e][0];
	    v1 = edgeArray[e][1];
	    status0 = ptList[v0].status;
	    status1 = ptList[v1].status;
	    
	    if ((status0 == status1))
	      {  
		if (status0==kPtBelowConeVertex)
		    
		  {
		    //or both pts are below cone vertex will not intersec
		    continue;
		  }
	      }
	    else if ((status0 == kPtBelowConeVertex) &&
		    (status1 == kPtOutsideConeBelowAxis))
	      {
		continue;
	      }
	    else if ((status1 == kPtBelowConeVertex) && 
		    (status0 == kPtOutsideConeBelowAxis))
	      {
		continue;
	      }
	    
	    E = ptList[v1].pt - ptList[v0].pt;
	    
	    delta = ptList[v0].pt - coneVertex;
	    deltaDotA = delta.dot(frustumAxis);
	    AdotE = frustumAxis.dot(E);//ptList[v1].deltaDotA - deltaDotA;
	    magESquared = E.mag2();
	    c2 = AdotE*AdotE - cosThetaSqr*magESquared;
	    if (c2 >= 0)
	      continue;
	    EdotDelta = E.dot(delta);
	    
	    c1 = AdotE*deltaDotA - cosThetaSqr*EdotDelta;
	    if (status0 == kPtBelowConeVertex)
	      {
		if ((c1 > -c2) ||
		   (c2*deltaDotA > c1*AdotE))
		  continue;
	      }
	    else if (status1 == kPtBelowConeVertex)
	      {
		if ((c1 < 0)||
		   (c2*deltaDotA > c1*AdotE))
		  continue;
	      }
	    else if (c1 < 0 || c1 > -c2)
	      continue;
	    c0 = deltaDotA*deltaDotA - cosThetaSqr*delta.mag2();
	    if (c1*c1 < c0*c2)
	      continue;
	    return true;
	  }
	//at this point if the bbox intersects any part of the axis of the infinite
	//cone then we must intersect this is the final weird set of cases
	// to check.
	BoundingBox newBBox(newBox[0],newBox[1]);
	if (newBBox.containsPoint(center1) || 
	   newBBox.containsPoint(center2) ||
	   newBBox.containsPoint(coneVertex) ||
	   newBBox.intersectsRay(coneVertex,frustumAxis,&dist))
	  return true;
      }
    return false;
  }

  // Stolen almost verbatim from Wild Magic library (www.magic-software.com)
  bool Box3::testIntersection(Ray3 const &ray,
			      double *distance,double inflateFactor) const
  {
    Box3 inflateBox;
    Box3 const *testBox;
    if (inflateFactor != 0)
    {
      inflateBox = createInflatedBox(*this,inflateFactor);
      testBox=&inflateBox;
    }
    else
      testBox=this;
    if (testBox->testIntersection(ray.getOrigin()))
      {
	*distance = 0.0;
	return true;
      }
    
    Vector3d extent = testBox->getDiagonal();
    extent *= 0.5;

    Vector3d kDiff(ray.getOrigin() - testBox->getCenter());

    double fT0 = 0.0;
    double fT1 = std::numeric_limits<double>::max();
    if (!FindIntersection(kDiff, ray.getDirection(), extent, 
			  fT0, fT1))
      return false;
    
    *distance = fT0;
    return true;
  }

  bool Box3::isDegenerate() const
  {
    for (int i=0; i<3; i++)
      {
	if ( (max()(i)-min()(i)) < ALMOST_ZERO)
	  return true;
      }

    return false;
  }
}
