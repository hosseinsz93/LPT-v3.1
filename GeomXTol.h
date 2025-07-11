#ifndef _GEOMX_TOLER_H_
#define _GEOMX_TOLER_H_

#include "GeomVector.h"

namespace GeomX
{
  // Calculate x to the Nth power.  N MUST be an integer.
  // IntPwrConst<N>::value(x)

  template < int N > struct IntPrwConstN
  { static double value(double pow)
    { return ((N & 1) ? pow : 1.0) * IntPrwConstN<N/2>::value(pow * pow); } };

  template <> struct IntPrwConstN < 0 >
  { static double value(double pow) { return 1.0; } };

  template < int N, bool > struct IntPrwConstB
  { static double value(double x) { return IntPrwConstN<N>::value(x); } };

  template < int N > struct IntPrwConstB < N, false >
  { static double value(double x) { return 1.0 / IntPrwConstN<-N>::value(x); }};

  template < int N > struct IntPwrConst
  { static double value(double x)
    { return IntPrwConstB < N, N >= 0 >::value(x); } };



  template < class T > struct tol_c
  {
      tol_c(T const & x, double tol)
        : _x(x)
        , _tol(tol)
      { }

      T const & _x;
      double _tol;
  };

  template < class T > tol_c<T> Tol(T const & x,
                                    double tol = IntPwrConst<-15>::value(10.0))
                                                   { return tol_c<T>(x, tol); }

  // Specify the tolerance as a power of 10.
  template < class T, int N /* = -15 */ > tol_c<T> Tol10(T const & x)
                           { return tol_c<T>(x, IntPwrConst<N>::value(10.0)); }



  template < class T, class U >
  inline bool operator<(tol_c<T> const & a, U const & b)
                                               { return (a._x + a._tol < b); }

  template < class T, class U >
  inline bool operator>=(tol_c<T> const & a, U const & b)
                                               { return ! (a < b); }

  template < class T, class U >
  inline bool operator>(tol_c<T> const & a, U const & b)
                                               { return (a._x - a._tol > b); }

  template < class T, class U >
  inline bool operator<=(tol_c<T> const & a, U const & b)
                                               { return ! (a > b); }

  template < class T, class U >
  inline bool operator==(tol_c<T> const & a, U const & b)
                                               { return ! (a != b); }

  template < class T, class U >
  inline bool operator!=(tol_c<T> const & a, U const & b)
                                               { return (a < b || a > b); }



  template < class T, class U >
  inline bool operator<(U const & a, tol_c<T> const & b)
                                               { return (a + b._tol < b._x); }

  template < class T, class U >
  inline bool operator>=(U const & a, tol_c<T> const & b)
                                               { return ! (a < b); }

  template < class T, class U >
  inline bool operator>(U const & a, tol_c<T> const & b)
                                               { return (a - b._tol > b._x); }

  template < class T, class U >
  inline bool operator<=(U const & a, tol_c<T> const & b)
                                               { return ! (a > b); }

  template < class T, class U >
  inline bool operator==(U const & a, tol_c<T> const & b)
                                               { return ! (a != b); }

  template < class T, class U >
  inline bool operator!=(U const & a, tol_c<T> const & b)
                                               { return (a < b || a > b); }


  template < >
  inline bool operator==(tol_c< Vector3d > const & a, Vector3d const & b)
  {
    return Tol(a._x(0), a._tol) == b(0) &&
           Tol(a._x(1), a._tol) == b(1) &&
           Tol(a._x(2), a._tol) == b(2);
  }
}

#endif
