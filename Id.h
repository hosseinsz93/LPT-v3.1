#ifndef _ID_H
#define _ID_H

// When using this class, you need to make sure that ni and nj are set first,
// otherwise the i,j,k conversions won't work.  IMPORTANT: This class will only
// work with one set of ni,nj values at a time.  To calculate for more than one,
// you need to change ni,nj between uses for different grids.

#include "GeomId.h"

class Id
{
  public:
    static void setDim(int _ni, int _nj)
    {
      ni = _ni;
      nj = _nj;
    }

    Id(Id const & _id)
      : id(_id.id)
    { }

    Id(Geom::Id const & _id)
      : id(_id.getValue())
    { }

    Id(int i, int j, int k)
    {
      id = Geom::Id(i + ni * (j + nj * k) + 1);
    }

    Id & operator=(Id const & _id)
    {
      if (this != &_id)
        id = _id.id;
      return *this;
    }

    void ijk(int & i, int & j, int & k) const
    {
      int tmp0 = id.getValue() - 1;
      int tmp1 = tmp0 / ni;
      i = tmp0 - ni * tmp1;
      k = tmp1 / nj;
      j = tmp1 - nj * k;
    }

    operator Geom::Id() const { return id; }

  private:
    static int ni, nj;
    Geom::Id id;
};

#endif
