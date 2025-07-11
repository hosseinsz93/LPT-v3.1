#ifndef _RAY3_TRACK_H
#define _RAY3_TRACK_H

#include <memory>

#include "GeomCoordSystem.h"
#include "GeomRay3.h"

#include "variables.h"
#include "UserCtx.h"

class Ray3Track : public Geom::Ray3
{
  public:
    Ray3Track()
      : Geom::Ray3()
      , cs()
      , faceId(-1)
      , pId(-1)
    { }
  
    Ray3Track(Vector3d const & origin,
              Vector3d const & direction,
              std::shared_ptr< UserCtx > & _user, int _pId)
      : Geom::Ray3(origin, direction)
      , faceId(-1)
      , user(_user)
      , pId(_pId)
    {
      cs.setOrthoNormalSet(origin, direction);
    }

    Ray3 * clone() const override
    {
      return const_cast<Ray3Track *>(this);
    }

    void setOrigin(Vector3d const &origin)
    {
      Geom::Ray3::setOrigin(origin);
      cs.setOrthoNormalSet(origin, getDirection());
    }

    void setDirection(Vector3d const &direction)
    {
      Geom::Ray3::setDirection(direction);
      cs.setOrthoNormalSet(getOrigin(), direction);
    }

    Geom::CoordSystem const & getCS() const
    {
      return cs;
    }

    void setFaceId(int i)
    {
      faceId = i;
    }

    int getFaceId() const
    {
      return faceId;
    }

    void setpId(int _pId)
    {
      pId = _pId;
    }

    int getpId() const
    {
      return pId;
    }


  private:
    Geom::CoordSystem cs;
    int faceId;

    // This is strictly for debugging, current particle id
    std::shared_ptr< UserCtx > user; 
    int pId;
};

#endif
