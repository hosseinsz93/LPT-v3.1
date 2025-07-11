#include "GeomXEdgeCalculations.h"


namespace GeomX
{
  template < int N, typename T >
  void edgePointVectors(VectorNT const & pnt,
                        VectorNT const & p0,     VectorNT const & p1,
                        VectorNT       & fromp0, VectorNT       & fromp1)
  {
    fromp0 = fromp1 = VectorNT(0.0);

    VectorNT v = VectorNT(pnt - p0).unit();
    VectorNT u = VectorNT(p1  - p0).unit();
    fromp0(0) = u.dot(v);
    fromp0(1) = VectorNT(v - fromp0(0) * u).length();

    v = VectorNT(pnt - p1).unit();
    u *= -1.0;
    fromp1(0) = u.dot(v);
    fromp1(1) = VectorNT(v - fromp1(0) * u).length();
  }

}
