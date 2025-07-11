#include "GeomXUtilities.h"

#include <iomanip>


namespace GeomX
{
  void debug()
  {
    std::cout << "debug" << std::endl;
  }


  // Find the angle from v0 to v1.  The result will be between 0 (inclusive)
  // and 360 (exclusive) following the right hand rule from v0 to v1.
  double angle2D(Vector3d v0, Vector3d v1)
  {
    v0.normalize();
    v1.normalize();

    double x = v0.dot(v1);
    if (x > 1.0) x = 1.0;
    else if (x < -1.0) x = -1.0;

    double y = v0(0)*v1(1) - v0(1)*v1(0);
    if (y > 1.0) y = 1.0;
    else if (y < -1.0) y = -1.0;

    if (std::abs(x) <= almostZero() && std::abs(y) <= almostZero())
      return 0.0;

    return atan2Deg(y, x);
  }


  // I had to copy the code here to get std::vector < bool > to work.
  void rotate(std::vector < bool > & a, int k)
  {
    if (k < 0 || k >= (int)a.size())
    {
      k %= (int)a.size();
      if (k < 0) k += (int)a.size();
    }
    if (k == 0) return;

    int c=0, v;
    for (v = 0; c < (int)a.size(); v++)
    {
      int t = v, tp = v + k;
      bool tmp = a[v];
      c++;
      while (tp != v)
      {
        a[t] = a[tp];
        t = tp;
        tp += k;
        if (tp >= (int)a.size()) tp -= (int)a.size();
        c++;
      }
      a[t] = tmp;
    }
  }


  // Take a number to an integer power.
  double pow(double d, int exp)
  {
    bool invert = (exp < 0);
    if (invert)
      exp *= -1;

    double pow = d;
    double result = 1.0;
    for ( ; exp; exp >>= 1, pow *= pow)
      if (exp & 1)
        result *= pow;

    return invert ? 1.0 / result : result;
  }


  std::ostream & operator<<(std::ostream & os, _Sci _s)
  {
    std::ios_base::fmtflags tmp(os.flags());
    os << std::scientific << std::left;
    std::streamsize prec = os.precision(_s._prec);
    os.width(_s._width);
    os << _s._d;
    os << std::setprecision(prec);
    os.flags(tmp);
    return os;
  }
}
