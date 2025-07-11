#include "Geomorient2d_sos.h"
#include "Geompredicates.h"
#include "GeomPoint2.h"

namespace Geom
{
#define SIGN(x) ((x) > 0. ? 1 : -1)
#define ORIENT1D(a,b) ((a) > (b) ? 1 : (a) < (b) ? -1 : 0)

  static int sortPointIds(Point2* points, unsigned n)
  {
    int sign = 1;

    for (unsigned int i = 0; i < n - 1; i++)
      for (unsigned int j = 0; j < n - 1 - i; j++)
        if ( points[j+1].getId() < points[j].getId() ) {
          Point2 tmp = points[j];

          points[j] = points[j+1];
          points[j+1] = tmp;
          sign = - sign;
        }
    return sign;
  }


  double orient2d_sos(Point2 const &point1,
                      Point2 const &point2,
                      Point2 const &point3,
                      double *orient2dresult)
  {
    double *A=  (double*) &point1.position()(0);
    double *B = (double*) &point2.position()(0);
    double *C = (double*) &point3.position()(0);

    double o=orient2d(A, B, C);

    if (orient2dresult)
      *orient2dresult = o;

    if (o != 0.)
      {
        //cout << " o!= 0. " << o << endl;
        return SIGN(o);
      }

    Point2 p[3];
    p[0] = point1;
    p[1] = point2;
    p[2] = point3;

    int const sign = sortPointIds(p, 3);

    /* epsilon^1/4 */
    o = ORIENT1D(p[1].position()(0), p[2].position()(0));
    //cout << "  /* epsilon^1/4 */ " << o << endl;
    if (o != 0.)
      return -SIGN(o)*sign;

    /* epsilon^1/2 */
    o = ORIENT1D(p[1].position()(1), p[2].position()(1));
    //cout << "  /* epsilon^1/2 */ " << o << endl;
    if (o != 0.)
      return SIGN(o)*sign;

    /* epsilon */
    o = ORIENT1D(p[0].position()(0), p[2].position()(0));
    //cout << "  /* epsilon */ " << o << endl;
    if (o != 0.)
      return SIGN(o)*sign;

    /* epsilon^3/2 */
    //cout << "     /* epsilon^3/2 */  " << sign <<  endl;
    return sign;

  }
}
