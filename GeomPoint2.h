#ifndef _GEOM_POINT2_H_
#define _GEOM_POINT2_H_

#include "GeomVector.h"

#include "GeomId.h"

#include "GeomCross2d.h"


namespace Geom
{
  class Point2
  {
    public:
      explicit Point2()
      {}

      Point2(Point2 const &rhs)
        : _position(rhs._position)
        , _id(rhs._id)
      {}

      explicit Point2(Point2 rhs,
                      Id const &id)
        : _position(rhs._position)
        , _id(id)
      {}

      explicit Point2(Vector2d const &c,
                      Id const &id = Id())
        : _position(c)
        , _id(id)
      {}

      explicit Point2(double x, double y, Id const &id = Id())
        : _id(id)
      {
        _position(0) = x;
        _position(1) = y;
      }

      void setId(Id const &id) { _id = id; }
      Id const &getId() const  { return _id; }

      Vector2d const &position() const { return _position; }
      Vector2d       &position()       { return _position; }

      Vector2d operator-(Point2 const &rhs) const
        { return Vector2d(_position - rhs._position); }

      Point2 & operator+=(Vector2d const &rhs)
      {
        _position += rhs;
        return *this;
      }

    private:
      Vector2d _position;
      Id _id;
  };
}

#endif // _MKGEOM_POINT2_H_
