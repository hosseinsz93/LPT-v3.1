#ifndef _GEOMX_UTILITIES_H_
#define _GEOMX_UTILITIES_H_


#include <cassert>
#include <cctype>
#include <cmath>
#include <vector>

// #include "boost/math/special_functions/fpclassify.hpp"
// #include "boost/math/special_functions/round.hpp"
// #include "boost/math/special_functions/trunc.hpp"

#include "GeomVector.h"


namespace GeomX
{
  void debug();

  // inline Vector2d toVector2D(Vector3d const & x)
  //   { return Vector2d(x(0), x(1)); }

  inline Vector2d const & toVector2D(Vector3d const & x)
    { return reinterpret_cast<Vector2d const &>(x); }

  inline Vector3d toVector3D(Vector2d const & x)
    { return Vector3d(x(0), x(1), 0.0); }

  inline Vector3d zeroZ(Vector3d const & x)
    { return Vector3d(x(0), x(1), 0.0); }

  template < class T > inline T sqr(T a) { return a*a; }

  template < class T > inline T sign(T const & a, T const & b)
    { return (b >= 0.0) ? std::abs(a) : -std::abs(a); }

  template < typename PeteExpression0, typename PeteExpression1 >
  double dist2D(PeteExpression0 const & a, PeteExpression1 const & b)
  {
    Vector2d va, vb;
    assign(va, a);
    assign(vb, b);
    return Vector2d(vb - va).length();
  }

  template < typename PeteExpression0, typename PeteExpression1 >
  double dist(PeteExpression0 const & a, PeteExpression1 const & b)
  {
    Vector3d va, vb;
    assign(va, a);
    assign(vb, b);
    return Vector3d(vb - va).length();
  }

  template < typename PeteExpression0, typename PeteExpression1 >
  double distsqr(PeteExpression0 const & a, PeteExpression1 const & b)
  {
    Vector3d va, vb;
    assign(va, a);
    assign(vb, b);
    return Vector3d(vb - va).mag2();
  }

  template < typename PeteExpression0, typename PeteExpression1 >
  Vector3d uv(PeteExpression0 const & a, PeteExpression1 const & b)
  {
    Vector3d va, vb;
    assign(va, a);
    assign(vb, b);
    return Vector3d(vb - va).unit();
  }


  inline void tolower(std::string & str)
  {
    for (std::string::iterator i=str.begin(); i!=str.end(); ++i)
      *i = std::tolower(*i);
  }



  inline double almostZero() { return 1.0e-15; }

  /*
    Constants to 500+ digits.
    pi
    3.1415926535897932384626433832795028841971693993751058209749445923078 \
    164062862089986280348253421170679821480865132823066470938446095505822 \
    317253594081284811174502841027019385211055596446229489549303819644288 \
    109756659334461284756482337867831652712019091456485669234603486104543 \
    266482133936072602491412737245870066063155881748815209209628292540917 \
    153643678925903600113305305488204665213841469519415116094330572703657 \
    595919530921861173819326117931051185480744623799627495673518857527248 \
    912279381830119491

    pi / 180
    1.7453292519943295769236907684886127134428718885417254560971914401710 \
    091146034494436822415696345094822123044925073790592483854692275281012 \
    398474218934047117319168245015010769561697553581238605305168788691271 \
    172087032963589602642490187704350918173343939698047594019224158946968 \
    481378963297818112495229298469927814479531045416008449560904606967176 \
    196468710514390888951836280826780369563245260844119508941294762613143 \
    108844183845478429899625621072806214155969235444237497596399365292916 \
    06237743435006638e-2

    180 / pi
    5.7295779513082320876798154814105170 332405472466564321549160243861202 \
    847148321552632440968995851110944186223381632864893281448264601248315 \
    036068267863411942122526388097467267926307988702893110767938261442638 \
    263158209610460487020506444259656841120171912057738566280431284962624 \
    203376187937297623870790340315980719624089522045186205459923396314841 \
    906966220115126609691801514787637366923164107126774038514690165499594 \
    192515711986479435210661624389035202306756177796757113315683506205731 \
    31336015650134889856e+1
  */

#ifndef M_PIl
  #define M_PIl 3.1415926535897932384626433832795029L  /* pi */
#endif

  inline double pi() { return M_PIl; }

  inline double degreeToRadian(double degree)
    { return degree * /* M_PIl / 180.0;*/
                 static_cast<double>(1.7453292519943295769236907684886127e-2L);
    }
  inline double radianToDegree(double radian)
    { return radian * /* 180.0 / M_PIl; */
                 static_cast<double>(5.7295779513082320876798154814105170e+1L);
    }

  // Force an angle in degrees to be between 0 and 360.
  inline double degree_000_360(double degree)
  {
    while (degree < 0.0)
      degree += 360.0;
    while (degree >= 360.0)
      degree -= 360.0;
    return degree;
  }

  // Force an angle in degrees to be between -180 and 180.
  inline double degree_N180_P180(double degree)
  {
    while (degree <= -180.0)
      degree += 360.0;
    while (degree > 180)
      degree -= 360.0;
    return degree;
  }

  // Force an angle in radians to be between 0 and 2*pi.
  inline double radian_000_2pi(double radian)
  {
    while (radian < 0.0)
      radian += 2.0 * pi();
    while (radian >= 2.0 * pi())
      radian -= 2.0 * pi();
    return radian;
  }

  // Force an angle in radians to be between -pi and pi.
  inline double radian_Npi_Ppi(double radian)
  {
    while (radian <= -pi())
      radian += 2.0 * pi();
    while (radian > pi())
      radian -= 2.0 * pi();
    return radian;
  }

  inline double sinDeg(double degree)
  {
    degree = degree_000_360(degree);

    if (degree == 0.0 || degree == 180.0)
      return 0.0;
    else if (degree == 90.0)
      return 1.0;
    else if (degree == 270.0)
      return -1.0;

    return std::sin(degreeToRadian(degree));
  }

  inline double cosDeg(double degree)
  {
    degree = degree_000_360(degree);

    if (degree == 90.0 || degree == 270.0)
      return 0.0;
    else if (degree == 0.0)
      return 1.0;
    else if (degree == 180.0)
      return -1.0;

    return std::cos(degreeToRadian(degree));
  }

  // First derivative.
  inline double dsinDeg(double degree)
  {
    static double const fac = degreeToRadian(1.0);
    return fac * cosDeg(degree);
  }

  inline double dcosDeg(double degree)
  {
    static double const fac = degreeToRadian(1.0);
    return -1 * fac * sinDeg(degree);
  }

  // Second derivative.
  inline double d2sinDeg(double degree)
  {
    static double const fac2 = sqr(degreeToRadian(1.0));
    return -1 * fac2 * sinDeg(degree);
  }

  inline double d2cosDeg(double degree)
  {
    static double const fac2 = sqr(degreeToRadian(1.0));
    return -1 * fac2 * cosDeg(degree);
  }



  // Returns angle >= 0 and <= 180.
  inline double acosDeg(double x)
  {
    double angle = std::acos(x);
    angle = radianToDegree(angle);
    return angle;
  }

  inline double atan2Deg(double y, double x)
  {
    double angle = std::atan2(y, x);
    angle = radianToDegree(angle);
    angle = degree_000_360(angle);
    return angle;
  }

  // Find the angle from v0 to v1.  The result will be between 0 (inclusive)
  // and 360 (exclusive) following the right hand rule from v0 to v1.
  double angle2D(Vector3d v0, Vector3d v1);

  inline
  double angle2D(Vector3d const & v01, Vector3d const & v00,
                 Vector3d const & v11, Vector3d const & v10)
  {
    return angle2D(Vector3d(v01 - v00), Vector3d(v11 - v10));
  }

  inline
  Vector3d rotate90AroundPlusZ(Vector3d const & c)
    { return Vector3d(-c(1), c(0), c(2)); }




  inline Vector2d convToPolar(Vector2d const & c)
    { return Vector2d(c.length(), atan2Deg(c(1), c(0))); }

  inline Vector3d convToPolar(Vector3d const & c)
  {
    return
          Vector3d(toVector2D(c).length(), atan2Deg(c(1), c(0)), c(2));
  }

  inline Vector2d convToRect(Vector2d const & c)
    { return Vector2d(cosDeg(c(1)) * c(0), sinDeg(c(1)) * c(0)); }

  inline Vector3d convToRect(Vector3d const & c)
    { return Vector3d(cosDeg(c(1)) * c(0), sinDeg(c(1)) * c(0), c(2)); }


  // Calculate the mod function using doubles.  This works with mod math
  // so you get the "correct" result for negative a.
  inline double fmod(double a, double b, double * returnN = 0)
  {
    double n = a / b;
    n = floor(n);
    if (returnN)
      *returnN = n;
    double mod = a - n * b;
    // When a is negative and very small such that the size of a < b's
    // smallest significant digits, b and mod are equal.
    if (mod == b)
      mod = 0.0;
    assert(mod >= 0.0 && mod < b);
    return mod;
  }


  // Calculate an alternate mod function where the value must be between the
  // min and max value.
  inline double fmod(double a, double min, double max)
  {
    a -= min;
    a = fmod(a, max-min);
    return a + min;
  }


  // inline double my_round(double num) { return (boost::math::round)(num); }
  // inline double my_trunc(double num) { return (boost::math::trunc)(num); }
  // inline double my_isnan(double num) { return (boost::math::isnan)(num); }
  // inline double my_isinf(double num) { return (boost::math::isinf)(num); }

  // Take a number to an integer power.
  double pow(double d, int exp);



  // Rotate array a of size n such that k is moved to the beginning.
  template < class T > void rotate(T * a, int n, int k)
  {
    if (k < 0 || k >= n)
    {
      k %= n;
      if (k < 0) k += n;
    }
    if (k == 0) return;

    int c=0, v;
    for (v = 0; c < n; v++)
    {
      int t = v, tp = v + k;
      T tmp = a[v];
      c++;
      while (tp != v)
      {
        a[t] = a[tp];
        t = tp;
        tp += k;
        if (tp >= n) tp -= n;
        c++;
      }
      a[t] = tmp;
    }
  }


  template < class T > void rotate(std::vector < T > & a, int k)
  {
    rotate(&*a.begin(), a.size(), k);
  }


  // I had to copy the code here to get std::vector < bool > to work.
  void rotate(std::vector < bool > & a, int k);


  class Interpolate
  {
    public:
      Interpolate(double ratio)
        : _ratio(ratio)
      { }

      Interpolate(double begin0, double value, double end0)
        { _ratio = (value - begin0) / (end0 - begin0); }

      template < class T >
      T operator()(T const & begin1, T const & end1)
        { return (T)(begin1 * (1.0 - _ratio) + end1 * _ratio); }

    private:
      double _ratio;
  };


  class Interpolate2
  {
    public:
      Interpolate2(double begin, double end)
        : _begin(begin)
        , _end(end)
      { }

      template < typename T >
      T operator()(double value, T const & begin, T const & end)
      {
        double ratio = (value - _begin) / (_end - _begin);
        return T(begin * (1.0 - ratio) + end * ratio);
      }

    private:
      double _begin;
      double _end;
  };


  // Given two ranges, find the secondary value given the primary value.
  class Ratio
  {
    public:
      Ratio(double priStart, double priEnd, double secStart, double secEnd)
        : _priStart(priStart)
        , _priEnd(priEnd)
        , _secStart(secStart)
        , _secEnd(secEnd)
      { }

      // Given the primary, interpolate the secondary.
      double operator()(double pri)
      {
        double ratio = (pri - _priStart) / (_priEnd - _priStart);
        return _secStart * (1.0 - ratio) + _secEnd * ratio;
      }

      // Given the secondary, interpolate the primary.
      double inverse(double sec)
      {
        double ratio = (sec - _secStart) / (_secEnd - _secStart);
        return _priStart * (1.0 - ratio) + _priEnd * ratio;
      }

    private:
      double _priStart;
      double _priEnd;
      double _secStart;
      double _secEnd;
  };


  class _Sci
  {
    public:
      inline _Sci(double d, int width, int prec)
        : _d(d)
        , _width(width)
        , _prec(prec)
      {}

      double _d;
      int _width;
      int _prec;
  };

  inline _Sci sci(double d, int width = 12, int prec = 6)
  {
    return _Sci(d, width, prec);
  }

  std::ostream & operator<<(std::ostream & os, _Sci _s);
}

#endif
