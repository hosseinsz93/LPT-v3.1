#ifndef _GEOM_LINE3
#define _GEOM_LINE3

#include <cmath>

#include "GeomVector.h"

namespace Geom
{
  class Segment3;

  /**
   * A line is the same as a Ray except for the fact that it is defined in both
   * directions (think of the intersection of 2 planes)
   */

  class Line3
  {
  public:
    Line3()
    {}

    Line3(Vector3d const &origin,
	  Vector3d const &direction) 
      : _origin(origin)
      , _direction(direction)
    {}

    Vector3d const & getOrigin() const {
      return _origin;
    }
    void setOrigin(Vector3d const &origin) {
      _origin = origin;
    }


    Vector3d const &getDirection() const {
      return _direction;
   }
    void setDirection(Vector3d const &direction) {
      _direction = direction;
    }
    double getDistanceSquared(Line3 const &line2,
			      double* pfThisLine=0, double* pfline2=0) const;

    /** Compute the distance between a line and a point. */
    double getDistanceSquared(Vector3d const &point,
			      double* linloc = 0) const;

    double getDistance(Vector3d const &point,
		       double* linloc = 0) const {
      return sqrt(getDistanceSquared(point,linloc));
    }

    /** Compute the shortest distance between a line and a segment. */
    double getDistanceSquared(Segment3 const &segment,
			      double* pfLinP = 0,
			      double* pfSegP = 0) const;
  
    double getDistance(Segment3 const &segment,
		       double* pfLinP = 0,
		       double* pfSegP = 0) const {
      return sqrt(getDistanceSquared(segment,pfLinP,pfSegP));
    }

  private:
    Vector3d _origin;
    Vector3d _direction;
  };
}

#endif
