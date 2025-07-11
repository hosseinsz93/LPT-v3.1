#ifndef _Geom_Vector_h_
#define _Geom_Vector_h_

#include <cmath>
#include <iostream>
#include <string>

#include "petsc.h"

#include "cmpnts.h"
#include "GeomCompileTimeAssert.h"


namespace Geom
{

template <int N, typename T, int I> struct VectorUnroll;

enum UninitializedToken {
  UNINITIALIZED,
};


template <int N, typename T> class Vector;

template <int N, typename T> Vector<N,T> operator*(double, Vector<N,T> const &);

template <int N, typename T>
class Vector
{
private:
  T v[N];

public:
  enum { SIZE = N };
  enum { NEL = N };

  /// const element access.
  template <typename I>
  inline T const &operator()(I const &i) const { return v[i]; }

  /// element access.
  template <typename I>
  inline T       &operator()(I const &i)       { return v[i]; }

  /// element access. operator() overloaded to permit indexing with enum type.
  inline T const &operator()(int i) const { return v[i]; }
  inline T       &operator()(int i)       { return v[i]; }

  /**
   * Default constructor, initializes to zero.
   */
  inline Vector()
    : v() // note, the empty parens cause zero initialization of primitive types.
  {}

  /**
   * Non-initializing constructor.
   */
  inline Vector(UninitializedToken const &)
  {}


  inline Vector(double const &x)
  { VectorUnroll<N,T,N-1>::assign(*this, x); }

  inline Vector(float const &x)
  { VectorUnroll<N,T,N-1>::assign(*this, x); }

  template <typename S>
  inline Vector(Vector<N,S> const &x)
  { VectorUnroll<N,T,N-1>::assign(*this, x); }

  /// 2D element-wise convenience constructor.
  template <typename T1,typename T2>
  inline Vector(T1 const &x0, T2 const &x1)
  { CompileTimeAssert(N == 2); (*this)(0) = x0; (*this)(1) = x1; }

  /// 3D element-wise convenience constructor.
  template <typename T1,typename T2,typename T3>
  inline Vector(T1 const &x0, T2 const &x1, T3 const &x2)
  { CompileTimeAssert(N == 3); (*this)(0) = x0; (*this)(1) = x1; (*this)(2) = x2; }

  /// 4D element-wise convenience constructor.
  template <typename T1,typename T2,typename T3,typename T4>
  inline Vector(T1 const &x0, T2 const &x1, T3 const &x2, T4 const &x3)
  { CompileTimeAssert(N == 4); (*this)(0) = x0; (*this)(1) = x1; (*this)(2) = x2; (*this)(3) = x3; }

  /// Convert Cmpnts.
  inline Vector(int i, int j, int k, Cmpnts const * const * const * const coor)
  {
    CompileTimeAssert(N == 3);
    Cmpnts const & c = coor[k][j][i];
    (*this)(0) = c.x;
    (*this)(1) = c.y;
    (*this)(2) = c.z;
  }
    
  inline Vector(Cmpnts const & coor)
  {
    CompileTimeAssert(N == 3);
    (*this)(0) = coor.x;
    (*this)(1) = coor.y;
    (*this)(2) = coor.z;
  }

  /// direct access to underlying array.
  T const (&asArray() const)[N] { return v; }
  T       (&asArray()      )[N] { return v; }

  /// == operator
  inline bool operator==(Vector<N, T> const &rhs) const {
    return equals(rhs);
  }
  inline bool operator!=(Vector<N, T> const &rhs) const {
    return !equals(rhs);
  }



  // Operators
  template <typename S>
  Vector & operator=(Vector<N,S> const & a)
  {
    if (this !=  (Vector<N,T> const *) &a)
      VectorUnroll<N,T,N-1>::assign(*this, a);
    return *this;
  }

  Vector operator+(Vector const & v) const
  {
    Vector r(UNINITIALIZED);
    VectorUnroll<N,T,N-1>::add(r, *this, v);
    return r;
  }

  Vector & operator+=(Vector const & x)
  {
    VectorUnroll<N,T,N-1>::addeq(*this, x);
    return *this;
  }

  Vector & operator+=(double const & x)
  {
    VectorUnroll<N,T,N-1>::addeq(*this, x);
    return *this;
  }

  Vector operator-() const
  {
    Vector r(UNINITIALIZED);
    VectorUnroll<N,T,N-1>::neg(r, *this);
    return r;
  }

  Vector operator-(Vector const & v) const
  {
    Vector r(UNINITIALIZED);
    VectorUnroll<N,T,N-1>::sub(r, *this, v);
    return r;
  }

  Vector & operator-=(Vector const & v)
  {
    VectorUnroll<N,T,N-1>::subeq(*this, v);
    return *this;
  }

  template < typename S >
  Vector operator*(S const & x) const
  {
    Vector r(UNINITIALIZED);
    VectorUnroll<N,T,N-1>::mult(r, *this, x);
    return r;
  }

  template < typename S >
  Vector const & operator*=(S const & x)
  {
    VectorUnroll<N,T,N-1>::multeq(*this, x);
    return *this;
  }

  template < typename S >
  Vector operator/(S const & x) const
  {
    Vector r(UNINITIALIZED);
    VectorUnroll<N,T,N-1>::div(r, *this, x);
    return r;
  }

  template < typename S >
  Vector const & operator/=(S const & x)
  {
    VectorUnroll<N,T,N-1>::diveq(*this, x);
    return *this;
  }



  /// dot product with another Vector.
  template <typename S>
  inline S dot(Vector<N,S> const &x) const
  { return VectorUnroll<N,T,N-1>::dot(*this,x); }

  /// Vector magnitude squared (sum of the square of the elements).
  inline T mag2() const { return dot(*this); }

  /// Vector magnitude (square root of mag2()).
  inline T mag() const { return sqrt(mag2()); }

  inline T length() const { return mag(); }

  /** Distance squared between this vector and another. */
  inline double distanceSquared(Vector<N,T> const &v) const {
    return Vector<N,T>(v - (*this)).mag2();
  }

  /** Distance between this vector and another. */
  inline double distance(Vector<N,T> const &v) const {
    return Vector<N,T>(v - (*this)).mag();
  }

  /// Unit Vector.
  inline Vector<N,T> unit() const
  {
    T const m = mag();
    return (m > T(0)) ? Vector<N,T>((*this)/m) : (*this);
  }

   /** Unitize vector and return scale factor. */
  inline T normalize() {
    T const len = mag(); // underflows to zero
    if (len == T(0)) {
      (*this) = T(0);
    }
    else {
      (*this) /= len;
    }
    return len;
  }

  /// Sum of the elements.
  inline T sum() const
  { return VectorUnroll<N,T,N-1>::sum(*this); }

  /// Max of the elements.
  inline T max() const
  { return VectorUnroll<N,T,N-1>::max(*this); }

  /// Min of the elements.
  inline T min() const
  { return VectorUnroll<N,T,N-1>::min(*this); }

  /// Element-wise maximum of this vector and another vector; output stored to result.
  inline void max(Vector<N,T> const &rhs, Vector<N,T> &result) const
  { VectorUnroll<N,T,N-1>::max(*this, rhs, result); }

  /// Element-wise minimum of this vector and another vector; output stored to result.
  inline void min(Vector<N,T> const &rhs, Vector<N,T> &result) const
  { VectorUnroll<N,T,N-1>::min(*this, rhs, result); }

  /// Product of elements
  inline T product() const
  { return VectorUnroll<N,T,N-1>::product(*this); }

  /// Is this numerically identical to another vector.
  inline bool equals(Vector<N,T> const &other) const
  { return VectorUnroll<N,T,N-1>::equals(*this, other); }

  /// P1-norm: sum of the absolute value of the elements
  inline T norm1() const
  { return VectorUnroll<N,T,N-1>::norm1(*this); }

  /// Average of the elements.
  inline T average() const { return sum()/N; }

  /// cross product with another Vector
  template <typename S>
  inline Vector<N,S> cross(Vector<N,S> const &x) const
  {
    CompileTimeAssert(N == 3);
    return Vector<N,S>
      ((*this)(1)*x(2) - (*this)(2)*x(1),
       (*this)(2)*x(0) - (*this)(0)*x(2),
       (*this)(0)*x(1) - (*this)(1)*x(0));
  }
};


template <int N, typename T>
std::ostream &
operator<<(std::ostream &output, Vector<N,T> const &x)
{
  unsigned int i;

  for (i = 0; i < 1; ++i)
    output << "(" << x(i);
  for (i = 1; i < N; ++i)
    output << " " << x(i);
  output << ")";
  return output;
}

/**
 * Template metaprogram to unroll Vector loops.
 */
template <int N, typename T, int I>
struct VectorUnroll
{
  static inline void
  assign(Vector<N,T> & v, float const & x)
  {
    VectorUnroll<N,T,I-1>::assign(v, x);
    v(I) = x;
  }

  static inline void
  assign(Vector<N,T> & v, double const & x)
  {
    VectorUnroll<N,T,I-1>::assign(v, x);
    v(I) = x;
  }

  template <typename S>
  static inline void
  assign(Vector<N,T> & v, Vector<N,S> const & x)
  {
    VectorUnroll<N,T,I-1>::assign(v, x);
    v(I) = x(I);
  }

  static inline void
  add(Vector<N,T> & r, Vector<N,T> const & a, Vector<N,T> const & b)
  {
    VectorUnroll<N,T,I-1>::add(r, a, b);
    r(I) = a(I) + b(I);
  }

  static inline void
  addeq(Vector<N,T> & v, Vector<N,T> const & x)
  {
    VectorUnroll<N,T,I-1>::addeq(v, x);
    v(I) += x(I);
  }

  static inline void
  addeq(Vector<N,T> & v, double const & x)
  {
    VectorUnroll<N,T,I-1>::addeq(v, x);
    v(I) += x;
  }

  static inline void
  neg(Vector<N,T> & r, Vector<N,T> const & a)
  {
    VectorUnroll<N,T,I-1>::neg(r, a);
    r(I) = -a(I);
  }

  static inline void
  sub(Vector<N,T> & r, Vector<N,T> const & a, Vector<N,T> const & b)
  {
    VectorUnroll<N,T,I-1>::sub(r, a, b);
    r(I) = a(I) - b(I);
  }

  static inline void
  subeq(Vector<N,T> & v, Vector<N,T> const & x)
  {
    VectorUnroll<N,T,I-1>::subeq(v, x);
    v(I) -= x(I);
  }

  template < typename S >
  static inline void
  div(Vector<N,T> & r, Vector<N,T> const & a, S const & b)
  {
    VectorUnroll<N,T,I-1>::div(r, a, b);
    r(I) = a(I) / b;
  }

  template < typename S >
  static inline void
  mult(Vector<N,T> & r, Vector<N,T> const & a, S const & b)
  {
    VectorUnroll<N,T,I-1>::mult(r, a, b);
    r(I) = a(I) * b;
  }

  template <typename S>
  static inline void
  multeq(Vector<N,T> & v, S & x)
  {
    VectorUnroll<N,T,I-1>::multeq(v, x);
    v(I) *= x;
  }

  template <typename S>
  static inline void
  diveq(Vector<N,T> & v, S & x)
  {
    VectorUnroll<N,T,I-1>::diveq(v, x);
    v(I) /= x;
  }





  template <typename S>
  static inline S
  dot(Vector<N,T> const &v1, Vector<N,S> const &v2)
  {
    return VectorUnroll<N,T,I-1>::dot(v1,v2) + v1(I)*v2(I);
  }

  static inline T
  sum(Vector<N,T> const &v)
  {
    return VectorUnroll<N,T,I-1>::sum(v) + v(I);
  }

  static inline T
  product(Vector<N,T> const &v)
  {
    return VectorUnroll<N,T,I-1>::product(v) * v(I);
  }

  static inline T
  max(Vector<N,T> const &v)
  {
    return std::max(VectorUnroll<N,T,I-1>::max(v), v(I));
  }

  static inline T
  min(Vector<N,T> const &v)
  {
    return std::min(VectorUnroll<N,T,I-1>::min(v), v(I));
  }

  static inline void
  min(Vector<N,T> const &v0, Vector<N,T> const &v1, Vector<N,T> &result)
  {
    result(I) = std::min(v0(I), v1(I));
    VectorUnroll<N,T,I-1>::min(v0, v1, result);
  }

  static inline void
  max(Vector<N,T> const &v0, Vector<N,T> const &v1, Vector<N,T> &result)
  {
    result(I) = std::max(v0(I), v1(I));
    VectorUnroll<N,T,I-1>::max(v0, v1, result);
  }

  static inline bool
  equals(Vector<N,T> const &v1, Vector<N,T> const &v2)
  {
    return VectorUnroll<N,T,I-1>::equals(v1,v2) && (v1(I) == v2(I));
  }

  static inline T
  norm1(Vector<N,T> const &v)
  {
    return VectorUnroll<N,T,I-1>::norm1(v) + abs(v(I));
  }
};


template <int N, typename T>
struct VectorUnroll<N,T,0>
{
  static inline void
  assign(Vector<N,T> & v, float const & x)
  {
    v(0) = x;
  }

  static inline void
  assign(Vector<N,T> & v, double const & x)
  {
    v(0) = x;
  }

  template <typename S>
  static inline void
  assign(Vector<N,T> & v, Vector<N,S> const & x)
  {
    v(0) = x(0);
  }

  static inline void
  add(Vector<N,T> & r, Vector<N,T> const & a, Vector<N,T> const & b)
  {
    r(0) = a(0) + b(0);
  }

  static inline void
  addeq(Vector<N,T> & v, Vector<N,T> const & x)
  {
    v(0) += x(0);
  }

  static inline void
  addeq(Vector<N,T> & v, double const & x)
  {
    v(0) += x;
  }

  template < typename S >
  static inline void
  div(Vector<N,T> & r, Vector<N,T> const & a, S const & b)
  {
    r(0) = a(0) / b;
  }

  template < typename S >
  static inline void
  mult(Vector<N,T> & r, Vector<N,T> const & a, S const & b)
  {
    r(0) = a(0) * b;
  }

  static inline void
  neg(Vector<N,T> & r, Vector<N,T> const & a)
  {
    r(0) = -a(0);
  }

  static inline void
  sub(Vector<N,T> & r, Vector<N,T> const & a, Vector<N,T> const & b)
  {
    r(0) = a(0) - b(0);
  }

  static inline void
  subeq(Vector<N,T> & v, Vector<N,T> const & x)
  {
    v(0) -= x(0);
  }

  template <typename S>
  static inline void
  multeq(Vector<N,T> & v, S & x)
  {
    v(0) *= x;
  }

  template <typename S>
  static inline void
  diveq(Vector<N,T> & v, S & x)
  {
    v(0) /= x;
  }





  template <typename S>
  static inline S
  dot(Vector<N,T> const &v1, Vector<N,S> const &v2)
  {
    return v1(0)*v2(0);
  }

  static inline T
  sum(Vector<N,T> const &v)
  {
    return v(0);
  }

  static inline T
  product(Vector<N,T> const &v)
  {
    return v(0);
  }

  static inline T
  max(Vector<N,T> const &v)
  {
    return v(0);
  }

  static inline T
  min(Vector<N,T> const &v)
  {
    return v(0);
  }

  static inline void
  min(Vector<N,T> const &v0, Vector<N,T> const &v1, Vector<N,T> &result)
  {
    result(0) = std::min(v0(0),v1(0));
  }

  static inline void
  max(Vector<N,T> const &v0, Vector<N,T> const &v1, Vector<N,T> &result)
  {
    result(0) = std::max(v0(0),v1(0));
  }

  static inline bool
  equals(Vector<N,T> const &v1, Vector<N,T> const &v2)
  {
    return v1(0) == v2(0);
  }

  static inline T
  norm1(Vector<N,T> const &v)
  {
    return abs(v(0));
  }
};

template <int N, typename T>
Vector<N,T> operator*(double s, Vector<N,T> const & v)
{
  Vector<N,T> r;
  VectorUnroll<N,T,N-1>::mult(r, v, s);
  return r;
}

}

typedef Geom::Vector<2,double> Vector2d;
typedef Geom::Vector<2,float>  Vector2f;
typedef Geom::Vector<3,double> Vector3d;
typedef Geom::Vector<3,float>  Vector3f;
// template <int N, typename T> typename Geom::Vector<N,T> VectorNT;
#define VectorNT Geom::Vector<N,T>

#endif
