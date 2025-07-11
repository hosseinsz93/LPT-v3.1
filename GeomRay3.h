#ifndef _GEOM_RAY3_H
#define _GEOM_RAY3_H

#include "GeomVector.h"

namespace Geom
{
  class Ray3
  {
  public:
    Ray3()
    {}

    Ray3(Vector3d const &origin,
	 Vector3d const &direction) 
      : _origin(origin)
      , _direction(direction)
    {}

    // I needed this route so I wouldn't loose the Ray3Track'ness if it was indeed a Ray3Track class.
    virtual Ray3 * clone() const
    {
      return const_cast<Ray3 *>(this);
    }

    Vector3d const & getOrigin() const
    {
      return _origin;
    }
    void setOrigin(Vector3d const &origin)
    {
      _origin = origin;
    }


    Vector3d const & getDirection() const
    {
      return _direction;
    }
    void setDirection(Vector3d const &direction)
    {
      _direction = direction;
    }

  private:
    Vector3d _origin;
    Vector3d _direction;
  };
}

#endif
