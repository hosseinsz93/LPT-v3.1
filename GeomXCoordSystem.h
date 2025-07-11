#ifndef _GEOMX_COORD_SYSTEM_H_
#define _GEOMX_COORD_SYSTEM_H_

#include <iostream>
#include <string>

#include "GeomVector.h"

#include "GeomCoordSystem.h"


namespace GeomX
{
  /**
   * Coordinate system class functions - definitions and transformations.
   */
  class CoordSystem : public Geom::CoordSystem
  {
  public:
    CoordSystem() : Geom::CoordSystem() { }

    CoordSystem (CoordSystem const & cs)
      : Geom::CoordSystem(static_cast<Geom::CoordSystem>(cs))
    { }

    CoordSystem (Geom::CoordSystem const & cs) : Geom::CoordSystem(cs) { }

    CoordSystem (Vector3d const & arg_origin,
                 Vector3d const & x_axis,
                 Vector3d const & y_axis,
                 Vector3d const & z_axis)

      : Geom::CoordSystem(arg_origin, x_axis, y_axis, z_axis)
    { }

    // A constructor that does the same thing as setOrthoNormalSet.
    CoordSystem (Vector3d const &arg_origin,
                 Vector3d const &normal)
      : Geom::CoordSystem(arg_origin, normal)
    { }

    CoordSystem & operator=(Geom::CoordSystem const & cs)
    {
      if (this == &cs)
        return *this;

      Geom::CoordSystem::operator=(cs);
      return dynamic_cast < CoordSystem & >(*this);
    }

    Vector3d convToProstarRotations() const;

    int printdebug(int detail) const;
  };


  struct prtProRot
  {
    prtProRot(CoordSystem const & cs) : _cs(cs) { }
    CoordSystem const & _cs;
  };

  std::ostream & operator<<(std::ostream & s, prtProRot const & ppr);


  struct outputLocalCommand
  {
      outputLocalCommand(CoordSystem const & cs, int id,
                         std::string const & type)
        : _cs(cs)
        , _id(id)
        , _type(type)
      { }

      outputLocalCommand(Geom::CoordSystem const & cs, int id,
                         std::string const & type)
        : _cs(cs)
        , _id(id)
        , _type(type)
      { }

      CoordSystem const _cs;
      int _id;
      std::string const _type;
  };

  std::ostream & operator<<(std::ostream & s, outputLocalCommand const & olc);

}

#endif
