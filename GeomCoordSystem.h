#ifndef _GEOM_COORDSYSTEM_H_
#define _GEOM_COORDSYSTEM_H_

#include <iostream>
#include <string>

#include "GeomVector.h"


namespace Geom
{
  /**
   * struct to define the rotation matrix.
   */
  struct Rot
  {
    Vector3d x[3];

    Vector3d       & operator[](int i)       { return x[i]; }
    Vector3d const & operator[](int i) const { return x[i]; }
  };


  /**
   * Coordinate system class functions - definitions and transformations.
   */
  class CoordSystem
  {
  public:
    CoordSystem()
    {
      origin = Vector3d(0, 0, 0);
      rot[0] = Vector3d(1, 0, 0);
      rot[1] = Vector3d(0, 1, 0);
      rot[2] = Vector3d(0, 0, 1);
    }

    CoordSystem (Vector3d const & arg_origin,
		 Vector3d const & x_axis,
		 Vector3d const & y_axis,
		 Vector3d const & z_axis)
    {
      origin = arg_origin;
      rot[0] = x_axis;
      rot[1] = y_axis;
      rot[2] = z_axis;
    }

    // A constructor that does the same thing as setOrthoNormalSet.
    CoordSystem (Vector3d const &arg_origin,
		 Vector3d const &normal)
    {
      setOrthoNormalSet(arg_origin, normal);
    }

    CoordSystem & operator=(CoordSystem const & cs)
    {
      if (this != &cs)
      {
        origin = cs.origin;
        rot[0] = cs.rot[0];
        rot[1] = cs.rot[1];
        rot[2] = cs.rot[2];
      }
      return *this;
    }

    // virtual ~CoordSystem() { }

    // Sets up a orthonormal coordinate system using the origin provided.  The
    // direction of normal becomes the z axis.  x and y are found such that
    // they are normal to z and of unit length.  The results are not unique.
    void setOrthoNormalSet (Vector3d const &arg_origin,
			    Vector3d const &normal)
    {
      origin = arg_origin;
      rot[0] = 0.0;
      rot[1] = 0.0;
      rot[2] = normal;
      calcOrthoNormalSet();
    }

    Vector3d rotInto (Vector3d const &crdIn) const
    {
      return Vector3d (crdIn.dot(rot[0]),
                       crdIn.dot(rot[1]),
                       crdIn.dot(rot[2]));
    }

    Vector3d transInto (Vector3d const &crdIn) const
    {
      Vector3d tmpCoord(crdIn - origin);
      return Vector3d (tmpCoord.dot(rot[0]),
                       tmpCoord.dot(rot[1]),
                       tmpCoord.dot(rot[2]));
    }

    Vector3d transInto2D (Vector3d const &crdIn) const
    {
      Vector3d wkCoord(crdIn - origin);
      return Vector3d(wkCoord.dot(rot[0]),
                      wkCoord.dot(rot[1]),
                      0);
    }

    Vector3d projOnto2D (Vector3d const &crdIn) const
    {
      return Vector3d(crdIn.dot(rot[0]),
                      crdIn.dot(rot[1]),
                      0);
    }

    double transIntoZ (Vector3d const &crdIn) const
    {
      Vector3d wkCoord(crdIn - origin);
      return wkCoord.dot(rot[2]);
    }

    Vector3d rotOutOf(Vector3d const &crdIn) const
    {
      return Vector3d(
        rot[0](0)*crdIn(0) + rot[1](0)*crdIn(1) + rot[2](0)*crdIn(2),
        rot[0](1)*crdIn(0) + rot[1](1)*crdIn(1) + rot[2](1)*crdIn(2),
        rot[0](2)*crdIn(0) + rot[1](2)*crdIn(1) + rot[2](2)*crdIn(2));
    }

    Vector3d transOutOf(Vector3d const &crdIn) const
    {
      return Vector3d(
      origin(0) + rot[0](0)*crdIn(0) + rot[1](0)*crdIn(1) + rot[2](0)*crdIn(2),
      origin(1) + rot[0](1)*crdIn(0) + rot[1](1)*crdIn(1) + rot[2](1)*crdIn(2),
      origin(2) + rot[0](2)*crdIn(0) + rot[1](2)*crdIn(1) + rot[2](2)*crdIn(2));
    }

    // Apply the lhs (*this) coordinate system to rhs coordinate system so
    // we get one system such that:
    //   local = result trans global
    //   local = lhs trans rhs trans global
    // give the same local coordinate.
    CoordSystem chain(CoordSystem const & rhs) const;

    void setOrigin(Vector3d const & _origin)
      { origin = _origin; }
    Vector3d const & getOrigin()   const { return origin; }
    Rot      const & getRot()      const { return rot;    }
    Vector3d const & getRot(int i) const { return rot[i]; }
    Vector3d const & getRotZ()     const { return rot[2]; }

    Vector3d       & getOrigin()         { return origin; }
    Vector3d       & getRot(int i)       { return rot[i]; }


  protected:
    void calcOrthoNormalSet();

    Vector3d origin;
    Rot rot;
  };
}

#endif
