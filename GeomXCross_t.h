#ifndef _GEOMX_CROSS_T_H_
#define _GEOMX_CROSS_T_H_

#include "GeomVector.h"

#include "GeomUtility.h"

#include "GeomCross2d.h"

namespace GeomX
{
  // The cross type is needed since 2d and 3d results have different return
  // types.  For 3d, we can use Vector<3>.  For 2d, we use the equivalent of
  // the 3d cross function where the zs are zero.  This means that, in the
  // resulting cross, only the z value is non-zero resulting cross product
  // and can be stored in a double.
  template < int N, typename T > class Cross_t;

  template < typename T >
  class Cross_t<2,T>
  {
    public:
      Cross_t() {}
      Cross_t(double x) : _x(x) { }

      T getX() { return _x; }

    private:
      T _x;
  };

  template < typename T >
  class Cross_t<3,T>
  {
    public:
      Cross_t() : _x() {}
      Cross_t(double x) : _x(x) { }
      explicit Cross_t(Geom::Vector<3,T> const & x) : _x(x) {}

      Geom::Vector<3,T> const & getX() { return _x; }

    private:
      Geom::Vector<3,T> _x;
  };



  template < typename T >
  inline T dot_t(Cross_t<2,T> v1, Cross_t<2,T> v2)
    { return v1.getX()*v2.getX(); }

  template < typename T >
  inline T dot_t(Cross_t<3,T> v1, Cross_t<3,T> v2)
    { return v1.getX().dot(v2.getX()); }



  template < typename T >
  inline Cross_t<2,T>
  cross_t(Geom::Vector<2,T> const v1, Geom::Vector<2,T> const v2)
    { return Cross_t<2,T>(Geom::cross2d(v1, v2)); }

  template < typename T >
  inline Cross_t<3,T>
  cross_t(Geom::Vector<3,T> const v1, Geom::Vector<3,T> const v2)
    { return Cross_t<3,T>(v1.cross(v2)); }



#if 0

  // The cross type is needed since 2d and 3d results have different
  // return types.  For 3d, we can use Vector<3>.  For 2d, we use the equivalent of the 3d cross function where the zs are zero.  This means that, in the resulting cross, only the z value is non-zero resulting cross product and can be stored in a double.
  template < int N, typename T > class Cross_t;

  template < typename T >
  class Cross_t <2,T>
  {
  public:
    inline Cross_t() { _x = 0.0; }
    inline Cross_t(T x) { _x = x; }
    inline Cross_t(Cross_t const &rhs) { _x = rhs._x; }
    inline Cross_t & operator=(Cross_t const &rhs)
    {
      _x = rhs._x;
      return *this;
    }

    inline operator T() { return _x; }

    inline Cross_t & operator*=(T const scale)
    {
      _x *= scale;
      return *this;
    }

    inline Cross_t operator*(T const scale) const
    {
      return Cross_t(_x * scale);
    }

    inline T normalize(T tol = Geom::ALMOST_ZERO)
    {
      T tmp = abs(_x);
      _x = (_x > tol) ? 1.0 : ((_x < -tol) ? -1.0 : 0.0);
      return tmp;
    }

  private:
    T _x;
  };



  template < int N, typename T1, typename T2 >
  struct Ops
  {
    inline static void multiplyAssign(T1 result, T2 scale)
    {
      Ops<N-1,T1,T2>::multiplyAssign(result, scale);
      result[N-1] *= scale;
    }
  };

  template < typename T1, typename T2 >
  struct Ops<0,T1,T2>
  {
    inline static void multiplyAssign(T1, T2) {}
  };



  template < typename T >
  class Cross_t<3,T> : public Geom::Vector<3,T>
  {
  public:
    inline Cross_t(Geom::Vector<3,T> const &x) : Geom::Vector<3,T>(x) {}

    inline Cross_t(T x, T y, T z)
    {
      (*this)[0] = x;
      (*this)[1] = y;
      (*this)[2] = z;
    }


    /*! Unary multiply (scale). (operator*=) */
    inline Cross_t & operator*=(T scale)
    {
      Ops<3, Cross_t, T >::multiplyAssign(*this, scale);
      return *this;
    }
  };



  template < int N, typename T >
  Cross_t<N,T>
  cross(VectorNT const v1, VectorNT const v2);

  template < typename T >
  inline Cross_t <2,T>
  cross(Vector <2,T> const v1, Vector <2,T> const v2)
  {
    return v1[0] * v2[1] - v1[1] * v2[0];
  }

  template < typename T >
  inline Cross_t<3,T>
  cross(Geom::Vector<3,T> const v1, Geom::Vector<3,T> const v2)
  {
    typedef Cross_t<3,T> Cross_t_t;
    return Cross_t_t(v1[1] * v2[2] - v1[2] * v2[1],
                     v1[2] * v2[0] - v1[0] * v2[2],
                     v1[0] * v2[1] - v1[1] * v2[0]);
  }



  /*! Compute the cross products (v1 X v2) X v3. */
  template < int N, typename T >
  VectorNT crossCross_t(VectorNT const &v1,
                        VectorNT const &v2,
                        VectorNT const &v3);

  template < typename T >
  inline Vector <2,T> crossCross_t(Vector <2,T> const &v1,
                                   Vector <2,T> const &v2,
                                   Vector <2,T> const &v3)
  {
    Cross_t <2,T> cp = cross(v1, v2);
    Vector <2,T> tmp;
    tmp[0] = -cp*v3[1];
    tmp[1] =  cp*v3[0];
    return tmp;
  }

  template < typename T >
  inline Geom::Vector<3,T> crossCross_t(Geom::Vector<3,T> const &v1,
                                        Geom::Vector<3,T> const &v2,
                                        Geom::Vector<3,T> const &v3)
  {
    return cross(cross(v1, v2), v3);
  }

#endif

}

#endif
