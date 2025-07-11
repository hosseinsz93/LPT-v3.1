#ifndef _CALC_PARTICLE_ACCELERATION_
#define _CALC_PARTICLE_ACCELERATION_

#include "variables.h"

#include "RungeKutta.h"

class RKparms : public Object
{
  public:
    UserCtx * u;
    Particle * p;

    int pId;

    PetscReal time;
    PetscReal time_o;
    PetscReal time_n;
};

RKY calcParticleAcceleration(Object & o, double h, double t, RKY const & yn);

#endif
