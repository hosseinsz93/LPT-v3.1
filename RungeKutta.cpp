#include "variables.h"

#include "GeomUtility.h"
#include "GeomXUtilities.h"
#include "ReysToCd.h"

#include "calcParticleAcceleration.h"
#include "RungeKutta.h"


constexpr int s = 6;

constexpr double a[][s-1] =
{
  {       0.0,            0.0,            0.0,             0.0,           0.0         },
  {     1.0/4.0,          0.0,            0.0,             0.0,           0.0         },
  {     3.0/32.0,       9.0/32.0,         0.0,             0.0,           0.0         },
  {  1932.0/2197.0, -7200.0/2197.0,  7296.0/2197.0,        0.0,           0.0         },
  {   439.0/216.0,       -8.0,       3680.0/513.0,    -845.0/4104.0,      0.0         },
  {    -8.0/27.0,         2.0,      -3544.0/2565.0,   1859.0/4104.0,  -11.0/40.0      }
};

constexpr double b[] =
{
       16.0/135.0,        0.0,       6656.0/12825.0, 28561.0/56430.0,  -9.0/50.0,      2.0/55.0
};

constexpr double bstar[] =
{
       25.0/216.0,        0.0,       1408.0/2565.0,   2197.0/4104.0,   -1.0/5.0,         0.0 
};

constexpr double c[] =
{
          0.0,          1.0/4.0,        3.0/8.0,        12.0/13.0,        1.0,         1.0/2.0
};


RKY operator*(double b, RKY const & a) { return RKY(a.getv() * b, a.getx() * b, a.getOutside()); }


RKY RungeKutta(Object & o, double h, double tn, RKY const & yn, Func f)
{
  RKparms & rkp = reinterpret_cast < RKparms & >(o);

  RKY k1 = f(o, h, tn,       yn           );
  if (k1.getOutside())
    return k1;
  RKY k2 = f(o, h, tn + h/2, yn + (h/2)*k1);
  if (k2.getOutside())
    return k2;
  RKY k3 = f(o, h, tn + h/2, yn + (h/2)*k2);
  if (k3.getOutside())
    return k3;
  RKY k4 = f(o, h, tn + h,   yn + h *   k3);
  if (k4.getOutside())
    return k4;
  return yn + (h/6) * (k1 + k2*2 + k3*2 + k4);
}



RKY RungeKuttaAdaptive(Object & o, double h, double tn, RKY const & yn, Func f)
{
  RKparms & rkp = reinterpret_cast < RKparms & >(o);

  if (rkp.pId == 0 && rkp.time >= 0.2)
    int debug = 1;

  RKY k0 = f(o, h, tn,          yn);
  if (k0.getOutside())
    return k0;
  RKY k1 = f(o, h, tn + c[1]*h, yn + (a[1][0]*k0) * h);
  if (k1.getOutside())
    return k1;
  RKY k2 = f(o, h, tn + c[2]*h, yn + (a[2][0]*k0 + a[2][1]*k1) * h);
  if (k2.getOutside())
    return k2;
  RKY k3 = f(o, h, tn + c[3]*h, yn + (a[3][0]*k0 + a[3][1]*k1 + a[3][2]*k2) * h);
  if (k3.getOutside())
    return k3;
  RKY k4 = f(o, h, tn + c[4]*h, yn + (a[4][0]*k0 + a[4][1]*k1 + a[4][2]*k2 + a[4][3]*k3) * h);
  if (k4.getOutside())
    return k4;
  RKY k5 = f(o, h, tn + c[5]*h, yn + (a[5][0]*k0 + a[5][1]*k1 + a[5][2]*k2 + a[5][3]*k3 + a[5][4]*k4) * h);
  if (k5.getOutside())
    return k5;

  RKY yn1     = yn + h * (b[0]*k0     + b[1]*k1     + b[2]*k2     + b[3]*k3     + b[4]*k4     + b[5]*k5);
  RKY yn1star = yn + h * (bstar[0]*k0 + bstar[1]*k1 + bstar[2]*k2 + bstar[3]*k3 + bstar[4]*k4 + bstar[5]*k5);

  RKY en1 = yn1 - yn1star;

  if (rkp.pId == 1)
    printf("pId: %d  time: %5f  er: %.5e\n", rkp.pId, rkp.time, en1.getv().mag() + en1.getx().mag());

  return yn1;
}
