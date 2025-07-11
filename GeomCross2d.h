#ifndef _GEOM_CROSS2D_H_
#define _GEOM_CROSS2D_H_

#include "GeomVector.h"


namespace Geom
{
  // Functions to calculate cross product with optimizations for 2D which
  // means that z is always zero.  When both operands are Vector2d,
  // a double is returned which is the z result of the cross.  The x and y
  // portions are zero since the original z's are zero.  When one operand is
  // double and the other is Vector2d, the double is used as the z
  // coordinate and x and y are assumed to be zero.

  inline double cross2d(Vector2d const &v0, Vector2d const &v1)
  {
    return v0(0)*v1(1) - v1(0)*v0(1);
  }

  inline Vector2d cross2d(double v0, Vector2d const &v1)
  {
    return Vector2d(-v0*v1(1), v0*v1(0));
  }

  inline Vector2d cross2d(Vector2d const &v0, double v1)
  {
    return Vector2d(v0(1)*v1, -v0(0)*v1);
  }
}

#endif // _MKGEOM_POINT2_H_
