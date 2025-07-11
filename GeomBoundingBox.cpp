#include "GeomBoundingBox.h"
#include "GeomUtility.h"

#include <limits>
#include <iostream>

namespace Geom
{

  void BoundingBox::reset()
  {
    double const maxValue = std::numeric_limits<double>::max();
    _min() = Vector3d(maxValue);
    _max() = Vector3d(-maxValue);
  }

  void BoundingBox::reset(Vector3d const &point)
  {
    _min() = _max() = point;
  }

  void BoundingBox::reset(Vector3d const &minc,
                          Vector3d const &maxc)
  {
    _min() = minc;
    _max() = maxc;
  }

  double BoundingBox::distanceSquared(Vector3d const &coord) const
  {
    double distancesquared = 0.0;
    double delta;

    if (coord(0) < _min(0))
      {
        delta = coord(0) - _min(0);
        distancesquared += delta*delta;
      }
    else if (coord(0) > _max(0))
      {
        delta = coord(0) - _max(0);
        distancesquared += delta*delta;
      }

    if (coord(1) < _min(1))
      {
        delta = coord(1) - _min(1);
        distancesquared += delta*delta;
      }
    else if (coord(1) > _max(1))
      {
        delta = coord(1) - _max(1);
        distancesquared += delta*delta;
      }

    if (coord(2) < _min(2))
      {
        delta = coord(2) - _min(2);
        distancesquared += delta*delta;
      }
    else if (coord(2) > _max(2))
      {
        delta = coord(2) - _max(2);
        distancesquared += delta*delta;
      }

    return distancesquared;
  }

  namespace {
    // return directional distance between boxes, or zero if overlapping
    inline double distance1D(double const& bb1min, double const& bb1max,
                             double const& bb2min, double const& bb2max)
    {
      if (bb1min > bb2max)
        return bb1min - bb2max;
      else if (bb2min > bb1max)
        return bb2min - bb1max;
      return 0.0;
    }
  }

  double BoundingBox::distanceSquared(BoundingBox const &bb) const
  {
    Vector3d delta(distance1D(_min(0), _max(0), bb._min(0), bb._max(0)),
                   distance1D(_min(1), _max(1), bb._min(1), bb._max(1)),
                   distance1D(_min(2), _max(2), bb._min(2), bb._max(2)));
    return delta.mag2();
  }


  bool BoundingBox::intersects(BoundingBox const &rhs) const
  {
    if (rhs._min(0) > _max(0) ||
        rhs._max(0) < _min(0))
      return false;
    if (rhs._min(1) > _max(1) ||
        rhs._max(1) < _min(1))
      return false;
    if (rhs._min(2) > _max(2) ||
        rhs._max(2) < _min(2))
      return false;
    return true;
  }

  bool BoundingBox::contains(BoundingBox const &rhs) const
  {
    return ( (_min(0) <= rhs._min(0)) &&
             (_max(0) >= rhs._max(0)) &&
             (_min(1) <= rhs._min(1)) &&
             (_max(1) >= rhs._max(1)) &&
             (_min(2) <= rhs._min(2)) &&
             (_max(2) >= rhs._max(2)) );
  }

  void BoundingBox::expand(Vector3d const &coord)
  {
    for (int i=0; i<3; i++)
      {
        _min(i) = std::min(_min(i), coord(i));
        _max(i) = std::max(_max(i), coord(i));
      }
  }

  void BoundingBox::expand(BoundingBox const &bb)
  {
    for (int i=0; i<3; i++)
      {
        _min(i) = std::min(_min(i), bb._min(i));
        _max(i) = std::max(_max(i), bb._max(i));
      }
  }

  void BoundingBox::inflate(double const size)
  {
    for (int i=0; i<3; i++)
      {
        _min(i) -= size;
        _max(i) += size;
      }
  }

  void BoundingBox::inflate(Vector3d const &size)
  {
    for (int i=0; i<3; i++)
      {
        _min(i) -= size(i);
        _max(i) += size(i);
      }
  }

  namespace {
    // compute min and max projected distance from a triangle to the origin
    // along the specified cartesian axis.
    void projectTriangle(int axis,
                         Vector3d const &akV0,
                         Vector3d const &akV1,
                         Vector3d const &akV2,
                         double &rfMin,
                         double &rfMax)
    {
      rfMin = akV0(axis);
      rfMax = rfMin;
      double fDot = akV1(axis);
      if (fDot < rfMin)
        rfMin = fDot;
      else if (fDot > rfMax)
        rfMax = fDot;

      fDot = akV2(axis);
      if (fDot < rfMin)
        rfMin = fDot;
      else if (fDot > rfMax)
        rfMax = fDot;
    }


    // compute min and max projected distance from a triangle to the origin
    // along the specified arbitrary vector
    void projectTriangle(Vector3d const &rkD,
                         Vector3d const &akV0,
                         Vector3d const &akV1,
                         Vector3d const &akV2,
                         double &rfMin,
                         double &rfMax)
    {
      rfMin = rkD.dot(akV0);
      rfMax = rfMin;
      double fDot = rkD.dot(akV1);
      if (fDot < rfMin)
        rfMin = fDot;
      else if (fDot > rfMax)
        rfMax = fDot;

      fDot = rkD.dot(akV2);
      if (fDot < rfMin)
        rfMin = fDot;
      else if (fDot > rfMax)
        rfMax = fDot;
    }


    void projectBox(Vector3d const &rkD,
                    BoundingBox const &rkBox,
                    double &rfMin,
                    double &rfMax)
    {
      double fDdC = rkD.dot(rkBox.center());
      Vector3d diag = rkBox.diagonal();
      double fR = 0.5*(diag(0)*fabs(rkD(0)) +
                       diag(1)*fabs(rkD(1)) +
                       diag(2)*fabs(rkD(2)));
      rfMin = fDdC - fR;
      rfMax = fDdC + fR;
    }
  } // unnamed namespace

  /*! Check if this box intersects with the specified triangle . */
  bool BoundingBox::intersectsTriangle(Vector3d const &v1,
                                       Vector3d const &v2,
                                       Vector3d const &v3,
                                       double tolerance) const
  {
    double fMin0, fMax0, fMin1, fMax1;

    // compute triangle normal
    Vector3d akE[3];
    akE[0] = v2 - v1;
    akE[1] = v3 - v1;
    Vector3d kD = akE[0].cross(akE[1]);
    kD.normalize();

    fMin0 = kD.dot(v1);
    fMax0 = fMin0;

    // test direction of triangle normal
    projectBox(kD, *this, fMin1, fMax1);
    if (fMax1 - fMin0 < -tolerance || fMax0 - fMin1 < -tolerance)
      return false;

    // test direction of box faces
    for (int i=0; i<3; i++)
      {
        projectTriangle(i, v1, v2, v3, fMin0, fMax0);
        fMin1 = _min(i);
        fMax1 = _max(i);
        if ( fMax1 - fMin0 < -tolerance || fMax0 - fMin1 < -tolerance )
          return false;
      }

    //  test direction of triangle-box edge cross products
    akE[2] = akE[1] - akE[0];
    for (int i0 = 0; i0 < 3; i0++)
      {
        // this loop is unrolled to exploit AABB
        // box x axis
        kD(0) = 0.0;
	kD(1) = akE[i0](2);
	kD(2) = -akE[i0](1);
        kD.normalize();
        projectTriangle(kD, v1, v2, v3, fMin0, fMax0);
        projectBox(kD, *this, fMin1, fMax1);
        if ( fMax1 - fMin0 < -tolerance || fMax0 - fMin1 < -tolerance )
          return false;

        // box y axis
        kD(0) = -akE[i0](2);
	kD(1) = 0.0;
	kD(2) = akE[i0](0);
        kD.normalize();
        projectTriangle(kD, v1, v2, v3, fMin0, fMax0);
        projectBox(kD, *this, fMin1, fMax1);
        if ( fMax1 - fMin0 < -tolerance || fMax0 - fMin1 < -tolerance )
          return false;

        // box z axis
        kD(0) = akE[i0](1);
	kD(1) = -akE[i0](0);
	kD(2) = 0.0;
        kD.normalize();
        projectTriangle(kD, v1, v2, v3, fMin0, fMax0);
        projectBox(kD, *this, fMin1, fMax1);
        if ( fMax1 - fMin0 < -tolerance || fMax0 - fMin1 < -tolerance )
          return false;
      }

    return true;
  }

  bool BoundingBox::intersectsSphere(Vector3d const &point,
                                     double radius) const
  {
    Vector3d center = this->center();
    Vector3d extent = this->diagonal();
    extent *= 0.5;

    double fAx = fabs(point(0) - center(0));
    double fAy = fabs(point(1) - center(1));
    double fAz = fabs(point(2) - center(2));
    double fDx = fAx - extent(0);
    double fDy = fAy - extent(1);
    double fDz = fAz - extent(2);

    if ( fAx <= extent(0) )
      {
        if ( fAy <= extent(1) )
          {
            if ( fAz <= extent(2) )
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
            if ( fAz <= extent(2) )
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
        if ( fAy <= extent(1) )
          {
            if ( fAz <= extent(2) )
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
            if ( fAz <= extent(2) )
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
  }

  namespace {
    // Stolen almost verbatim from Wild Magic library (www.magic-software.com)
    bool clip (double fDenom, double fNumer, double &rfT0, double &rfT1)
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
    bool findIntersection(Vector3d const &point,
                          Vector3d const &direction,
                          Vector3d const &extent,
                          double &rfT0,
                          double &rfT1)
    {
      double fSaveT0 = rfT0, fSaveT1 = rfT1;

      bool bNotEntirelyClipped =
        clip(+direction(0),-point(0)-extent(0),rfT0,rfT1) &&
        clip(-direction(0),+point(0)-extent(0),rfT0,rfT1) &&
        clip(+direction(1),-point(1)-extent(1),rfT0,rfT1) &&
        clip(-direction(1),+point(1)-extent(1),rfT0,rfT1) &&
        clip(+direction(2),-point(2)-extent(2),rfT0,rfT1) &&
        clip(-direction(2),+point(2)-extent(2),rfT0,rfT1);

      return bNotEntirelyClipped && ( rfT0 != fSaveT0 || rfT1 != fSaveT1 );
    }
  } // unnamed namespace

  // Stolen almost verbatim from Wild Magic library (www.magic-software.com)
  bool BoundingBox::intersectsRay(Vector3d const &point,
                                  Vector3d const &direction,
                                  double *distance) const
  {
    if (this->containsPoint(point))
      {
        *distance = 0.0;
        return true;
      }

    Vector3d center = this->center();
    Vector3d extent = this->diagonal();
    extent *= 0.5;

    Vector3d kDiff(point - center);

    double fT0 = 0.0;
    double fT1 = std::numeric_limits<double>::max();
    if (!findIntersection(kDiff, direction, extent, fT0, fT1))
      return false;

    *distance = fT0;
    return true;
  }

  bool BoundingBox::intersectsSegment(Geom::Vector<3, double> const &p1,
                                      Geom::Vector<3, double> const &p2)
  { // define the ray direction
	if(this->containsPoint(p1) || this->containsPoint(p2))
	  return true;

        Geom::Vector<3, double> dir(p2(0)-p1(0), p2(1)-p1(1), p2(2)-p1(2) );
	double segDist = dir.mag();
    dir.normalize();
    double rayDist;
    if(!intersectsRay(p1, dir, &rayDist))
    	return false;

    if(rayDist < segDist)
    	return true;

    return false;
  }

  std::ostream & operator<< (std::ostream & out, BoundingBox const &bbox)
  {
    out << "("   << bbox.min() << " , " << bbox.max() << ")" << std::endl;
    return out;
  }
}
