#ifndef _RUNGEKUTTA_H_
#define _RUNGEKUTTA_H_

#include "variables.h"


class RKY
{
  public:
    RKY() { }
    RKY(Vector3d const & _v, Vector3d const & _x, bool _outside = 0) : v(_v), x(_x), outside(_outside) {}
    RKY const & operator=(RKY const & _y)
    {
      if (this != &_y)
      {
        v = _y.v;
        x = _y.x;
        outside = _y.outside;
      }
      return *this;
    }
    
    Vector3d const & getv() const { return v; }
    Vector3d const & getx() const { return x; }
    bool       getOutside() const { return outside; }
    RKY operator+(RKY const & b) const { return RKY(this->v + b.v, this->x + b.x, this->outside); }
    RKY operator-(RKY const & b) const { return RKY(this->v - b.v, this->x - b.x, this->outside); }
    RKY operator*(PetscScalar b) const { return RKY(this->v * b, this->x * b, this->outside); }
    RKY operator/(PetscScalar b) const { return RKY(this->v / b, this->x / b, this->outside); }

    double mag() const { return v.mag() + x.mag(); }
  
  private:
    // Particle velocity
    Vector3d v;

    // Particle location
    Vector3d x;

    // Is the partile outside the bakground mesh domain?
    int outside;
};

RKY operator*(double b, RKY const & a);

typedef RKY (Func)(Object & o, double h, double t, RKY const & yn);

// Integrate one timestep with Runge-Kutta.
RKY RungeKutta(Object & o, double h, double tn, RKY const & yn, Func f);

RKY RungeKuttaAdaptive(Object & _o, double h, double tn, RKY const & yn, Func f);

#endif  // _RUNGEKUTTA_H
