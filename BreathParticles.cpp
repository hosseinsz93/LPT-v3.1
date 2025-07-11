#include "BreathParticles.h"

#include <cmath>
#include <cstdlib>

#include "GeomXTol.h"


void breathAddParticles(double time, std::shared_ptr< UserCtx > & user,
                        std::vector<Particle> & Ptc)
{
  constexpr double ymin = 0.98;
  constexpr double ymax = 1.02;
  constexpr double ydel = ymax - ymin;
  constexpr double zmin = 1.67;
  constexpr double zmax = 1.68;
  constexpr double zdel = zmax - zmin;
  
  // Do we need to add particles?
  if (fmod(time, 5.0) == GeomX::Tol(0.0))
  {
    for (int i=0; i<6; ++i)
    {
      double r = (double)rand() / RAND_MAX;
      double y = ymin + r * ydel;

             r = (double)rand() / RAND_MAX;
      double z = zmin + r * zdel;

      Ptc.push_back(Particle());
      auto & p = Ptc.back();
      p.coor.x = 0.005;
      p.coor.y = y;
      p.coor.z = z;
      // p.ztraj.push_back(p.coor);
      p.stime = time;
      p.DDenPt = 997.0;  // Density of water kg/m^3
    }
    user->ptc_num = Ptc.size();
  }
}
  

// Calculate particle diameter.
void breathCalcPartDia(double time, Particle & P)
{
  P.etime = time;
  double t = time - P.stime;

#if 1
  // Calculation given to me on 2021-03-18 by Ali
  P.DDiaPt = 10.0e-6;
  if (t >= 0.0)
  {
    P.DDiaPt = ((-700e-6*t - 70e-6)*t - 5e-6)*t + 10e-6;
    if (P.DDiaPt < 0.1e-6)
      P.DDiaPt = 0.1e-6;
    if (P.DDiaPt > 10.0e-6)
      P.DDiaPt = 10.0e-6;
  }
#else
  // Original calculation
  if (t >= 0.0)
  {
    P.dia = (0.0011*t - 0.2133) * t + 10.1;
    P.dia /= 1000000.0;
  }
#endif
}

void breathCalcPartDia(double time, std::vector<Particle> & Ptc)
{
  for (auto & P : Ptc)
    breathCalcPartDia(time, P);
}
