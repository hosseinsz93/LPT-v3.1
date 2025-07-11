#include "calcParticleAcceleration.h"

#include "GeomUtility.h"
#include "GeomXUtilities.h"
#include "ReysToCd.h"


// Subroutines used from other parts of the program, mostly, if not all,
// in tracking.cpp
PetscErrorCode ParticleLocation(UserCtx * user, Particle *P,
                                PetscInt iradial, int pId, bool outsideIsAnError = true);

PetscTruth ISPointInCell(Cmpnts p, UserCtx * user,
                         Location loc, Cmpnts *coef, Particle *P, int pId,
                         bool printErr = true);

PetscErrorCode VelocityInterpolation(UserCtx * user,
                                     Particle *P, PetscReal time,
                                     PetscReal time_o, PetscReal time_n, bool VF0update = true);


RKY calcParticleAcceleration(Object & o, double h, double t, RKY const & y)
{
  RKparms & rkp = reinterpret_cast < RKparms & >(o);
  UserCtx * u = rkp.u;
  Particle * p = rkp.p;

  double const L = u->lengScale;
  double const U = u->velScale;
  double const T = L / U;

  // Get the fluid velocity at the current location of the particle.
  p->geomLookupStatus = GEOM_NOT_DONE;
  p->coor = Cmpnts(y.getx()(0), y.getx()(1), y.getx()(2));
  if (ParticleLocation(u, p, 5, rkp.pId/*, false*/) == OUTSIDE_UNKNOWN)
    return RKY(y.getv(), Vector3d(p->coor), true);

  Location sloc;
  ISPointInCell(p->coor, u, sloc, &p->coef, p, rkp.pId);

  VelocityInterpolation(u, p, rkp.time, rkp.time_o, rkp.time_n, false);


  // Calculate the drag force.
  Vector3d Vpf(y.getv() - p->Vf1);

  double rey = p->Rey = u->DDenFd * p->DDiaPt * (U*Vpf.mag()) / u->DVisFd;

  Vector3d drag;

  // The following true path stops divide by zero and gets a valid result.
  if (rey <= 0.001)
    drag = -3.0 * Geom::PI * u->DVisFd * p->DDiaPt * (U*Vpf);

  else
  {
    double Cd = ReysToCd::conv(rey);
    drag = -0.125 * u->DDenFd * Cd * Geom::PI
      * GeomX::pow(p->DDiaPt, 2) * (U*Vpf.mag()) * (U*Vpf);
  }

  // Calculate the particle mass.
  double DVolPt = Geom::PI * GeomX::pow(p->DDiaPt, 3) / 6;
  double DMasPt = p->DDenPt * DVolPt;

  // Calculate the total acceleration including gravity.
  Vector3d acc = drag / DMasPt + u->DAccGrav;

  // Calculate the terminal velocity at this acceleration.
  double termVel = acc.mag() * GeomX::pow(p->DDiaPt, 2) * (p->DDenPt - u->DDenFd)
                 / (18.0 * u->DVisFd);

  // If the acceleration converted to velocity is greater than the terminal velocity,
  // use the terminal velocity for acceleration (after an appropriate scaling).
  if (acc.mag() * h * T > termVel)
    acc = acc.normalize() * (termVel / (h * T));

  // non-dimensionalize acceleration.
  acc *= T/U;

  return RKY(y.getv(), acc);
}
