#include <math.h>
#include <stdio.h>

#include "GeomUtility.h"
#include "ReysToCd.h"


int main(int arg , char **argv)
{
  double dia = 30.0e-6;
  double den = 1.176;
  double vis = 1.86E-05;
  double vol = 4.0 / 3.0 * Geom::PI * pow(dia, 3);
  double mas = den * vol ;

  constexpr double powers[] = { 1e-6, 1e-5, 1e-4, 1e-3, 1e-2,
                                1e-1, 1e-0, 1e+1, 1e+2 };
  int pcnt = sizeof(powers) / sizeof(double);

  for (int i = 0; i < pcnt; ++i)
  {
    for (double rey=1.0; rey<10.0; rey+=1.0)
    {
      double Cd = 0.0;
    
      if (rey*powers[i] <= 5.0e-6)
      {
        Cd = 0.0;
      }

      else if (rey*powers[i] <= 0.001)
      {
        // DFdrag1 = -3.0 * Geom::PI * user->DVisFd * P->DDiaPt * (U*Vpf0);
        Cd = 3.0 * Geom::PI * vis * dia * 1.0 / mas;
      }

      else
      {
        Cd = ReysToCd::conv(rey*powers[i]);
      }


      double Cd2 = 0.0;
      if (rey*powers[i] > 1000)
        Cd2 = 0.424;
      else
        Cd2 = 24 / rey*powers[i] * (1. + 1./6. * pow(rey*powers[i], 2.0/3.0));

      
      printf("%e\t%e\t%e\n", rey*powers[i], Cd, Cd2);
    }
  }
}
