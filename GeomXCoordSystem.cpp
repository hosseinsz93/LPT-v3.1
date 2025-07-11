#include "GeomXCoordSystem.h"

#include <limits>

#include "GeomUtility.h"


namespace GeomX
{
/*
      Convert rotation matrix to axes rotations needed
      to define the rotation matrix in prostar.

Prostar defines the rotation matrix as a series of rotations around coordinate
axes in the following order:

  XY - right hand rotation around z axis - rotate x axis towards the y axis.
  YZ - right hand rotation around x axis - rotate y axis towards the z axis.
  ZX - right hand rotation around y axis - rotate z axis towards the x axis.

NOTE: Prostar stores each axis unit vector in column order.  There routines
      store them in row order.  You will need to transpose these matrices to
      get the number to match in prostar.

The rotation matrices are:

                         [  cos(XY)    sin(XY)   0.0   ]
                rotxy := [ -sin(XY)    cos(XY)   0.0   ]
                         [  0.0        0.0       1.0   ]

                         [  1.0        0.0       0.0   ]
                rotyz := [  0.0      cos(YZ)   sin(YZ) ]
                         [  0.0     -sin(YZ)   cos(YZ) ]

                         [  cos(ZX)    0.0    -sin(ZX) ]
                rotzx := [  0.0        1.0       0.0   ]
                         [  sin(ZX)    0.0     cos(ZX) ]

Prostar multiplies these matrices together this way:

  rotFull = rotzx (rotyz rotxy)

    to get the full rotation matrix.

rotFull :=
  [cos(XY) cos(ZX) - sin(ZX) sin(YZ) sin(XY) ,
                     sin(XY) cos(ZX) + sin(ZX) sin(YZ) cos(XY) ,
                                                              -sin(ZX) cos(YZ)]
  [-cos(YZ) sin(XY) ,cos(YZ) cos(XY) ,                         sin(YZ)        ]
  [cos(XY) sin(ZX) + cos(ZX) sin(YZ) sin(XY) ,
                     sin(XY) sin(ZX) - cos(ZX) sin(YZ) cos(XY) ,
                                                               cos(ZX) cos(YZ)]


Here's how to get the prostar rotations out of this matrix:

From the matrix, we have these equations:

  -sin(ZX) cos(YZ) = rotFull(0)(2)
  -cos(YZ) sin(XY) = rotFull(1)(0)
   cos(YZ) cos(XY) = rotFull(1)(1)
   cos(ZX) cos(YZ) = rotFull(2)(2)

After some algebra:

  sin(XY) = -rotFull(1)(0) / cos(YZ)
  cos(XY) =  rotFull(1)(1) / cos(YZ)

  sin(ZX) = -rotFull(0)(2) / cos(YZ)
  cos(ZX) =  rotFull(2)(2) / cos(YZ)



If cos(YZ) != 0, we can calculate XY and ZX as follows:

  XY = atan2(-rotFull(1)(0) / cos(YZ), rotFull(1)(1) / cos(YZ))
  ZX = atan2(-rotFull(0)(2) / cos(YZ), rotFull(2)(2) / cos(YZ))

Since cos(YZ) is in both terms, it doesn't matter to the result of the atan2
function giving this result:

  XY = atan2(-rotFull(1)(0), rotFull(1)(1))
  ZX = atan2(-rotFull(0)(2), rotFull(2)(2))

Now calculate the value of YZ:

  sin(XY) = -rotFull(1)(0) / cos(YZ)
  cos(XY) =  rotFull(1)(1) / cos(YZ)

becomes:

  cos(YZ) = -rotFull(1)(0) / sin(XY)
  cos(YZ) =  rotFull(1)(1) / cos(XY)

sin(XY) and cos(XY) are now known.  However, it is possible for either one,
but not both, values to be zero.  Use the equation with the largest value to
calculate YZ.

  YZ = atan2(sin(YZ), cos(YZ))
  YZ = atan2(rotFull(1)(2), cos(YZ))



What if cos(YZ) == 0?

This implies that sin(YZ) is +-1 or that YZ is +-90 degrees.  Substituting
these values into rotFull gives:

rotFull :=
  [cos(XY) cos(ZX)-+sin(ZX) sin(XY) , sin(XY) cos(ZX) +- sin(ZX) cos(XY) , 0.0]
  [0.0 ,                              0.0 ,                              +-1.0]
  [cos(XY) sin(ZX)+-cos(ZX) sin(XY) , sin(XY) sin(ZX) -+ cos(ZX) cos(XY) , 0.0]

Take the following equations out of rotFull:

  cos(XY) cos(ZX) -+ sin(ZX) sin(XY) = rotFull(0)(0)
  sin(XY) cos(ZX) +- sin(ZX) cos(XY) = rotFull(0)(1)

Apply the trigonometric identities:

  [[sine.sum.gif]]

  [[cosine.sum.gif]]

gives:

  sin(XY) cos(ZX) +- sin(ZX) cos(XY) = sin(XY +- ZX) = rotFull(0)(1)
  cos(XY) cos(ZX) -+ sin(XY) sin(ZX) = cos(XY +- ZX) = rotFull(0)(0)

  XY +- ZX = atan2(sin(XY +- ZX), cos(XY +- ZX))
  XY +- ZX = atan2(rotFull(0)(1), rotFull(0)(0))

Since XY and ZX are dependent and to ignore the sign of the sine, we are
going to set ZX to zero leaving:

  XY = atan2(rotFull(0)(1), rotFull(0)(0))

 */

  // Getting rid of negarive zero.  Just another example of compulsive
  // obsessive disorder.
  double removeNegZero(double x)
  {
    if (abs(x) < std::numeric_limits<double>::min())
      x = 0.0;
    return x;
  }

  Vector3d CoordSystem::convToProstarRotations() const
  {
    double xy, yz, zx;

    if (std::abs(rot[1](2)) < 0.999999)
    {
      xy = removeNegZero(
                   Geom::RADIAN_TO_DEGREE * std::atan2(-rot[1](0), rot[1](1)));
      zx = removeNegZero(
                   Geom::RADIAN_TO_DEGREE * std::atan2(-rot[0](2), rot[2](2)));

      double sinXY = std::sin(Geom::DEGREE_TO_RADIAN * xy);
      double cosXY = std::cos(Geom::DEGREE_TO_RADIAN * xy);
      double cosYZ;
      if (std::abs(sinXY) > std::abs(cosXY))
        cosYZ = -rot[1](0) / sinXY;
      else
        cosYZ =  rot[1](1) / cosXY;
      yz = removeNegZero(
                        Geom::RADIAN_TO_DEGREE * std::atan2(rot[1](2), cosYZ));
    }
    else
    {
      xy = removeNegZero(
                    Geom::RADIAN_TO_DEGREE * std::atan2(rot[0](1), rot[0](0)));
      yz = (rot[1](2) > 0.0) ? 90 : -90;
      zx = 0.0;
    }

    return Vector3d(xy, yz, zx);
  }


  int CoordSystem::printdebug(int detail) const
  {
    std::cout << "<< class CoordSystem >>"  << std::endl;
    std::cout << "origin: " << origin       << std::endl;
    std::cout << "rot x:  " << rot[0]       << std::endl;
    std::cout << "rot y:  " << rot[1]       << std::endl;
    std::cout << "rot z:  " << rot[2]       << std::endl;
    std::cout << "<< end of CoordSystem >>" << std::endl;
    return 0;
  }


  std::ostream & operator<<(std::ostream & s, prtProRot const & ppr)
  {
    Vector3d proRot = ppr._cs.convToProstarRotations();
    return s << proRot(0) << "," << proRot(1) << "," << proRot(2);
  }


  std::ostream & operator<<(std::ostream & s, outputLocalCommand const & olc)
  {
    return s << "local," << olc._id << "," << olc._type
             << ", " << olc._cs.getOrigin()(0)
             << ","  << olc._cs.getOrigin()(1) << "," << olc._cs.getOrigin()(2)
             << ", " << prtProRot(olc._cs) << std::endl;
  }

}
