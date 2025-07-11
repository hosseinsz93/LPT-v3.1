#ifndef _CMPNTS_h_
#define _CMPNTS_h_

class Cmpnts
{
  public:
    PetscScalar x, y, z;

    Cmpnts() : x(0.0), y(0.0), z(0.0) { }

    Cmpnts(PetscScalar _x, PetscScalar _y, PetscScalar _z) : x(_x), y(_y), z(_z) { }

    Cmpnts operator+(Cmpnts const & b) const
    {
      return Cmpnts(this->x + b.x,
                    this->y + b.y,
                    this->z + b.z);
    }

    Cmpnts operator*(PetscScalar mult) const
    {
      return Cmpnts (this->x * mult,
                     this->y * mult,
                     this->z * mult);
    }
};

#endif
