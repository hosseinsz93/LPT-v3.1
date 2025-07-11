#include "Geomorient3d_sos.h"
#include "Geompredicates.h"
#include "GeomPoint3.h"

namespace Geom
{
#define SIGN(x) ((x) > 0. ? 1 : -1)
#define ORIENT1D(a,b) ((a) > (b) ? 1 : (a) < (b) ? -1 : 0)

  static int sortPointIds(Point3* points, unsigned n)
  {
    int sign = 1;
    unsigned i, j;
    
    for (i = 0; i < n - 1; i++)
      for (j = 0; j < n - 1 - i; j++)
        if ( points[j+1].getId() < points[j].getId() ) {
          Point3 tmp = points[j];
          
          points[j] = points[j+1];
          points[j+1] = tmp;
          sign = - sign;
        }
    return sign;
  }





  double orient3d_sos(Point3 const &point1,
                      Point3 const &point2,
                      Point3 const &point3,
                      Point3 const &point4,
                      double *orient3dresult)
  {
    double *A=  (double*) &point1.position()(0);
    double *B = (double*) &point2.position()(0);
    double *C = (double*) &point3.position()(0);
    double *D = (double*) &point4.position()(0);

    double o=orient3d (A,B,C,D);

    if (orient3dresult)
      *orient3dresult = o;

    if (o != 0.)
      {
        //cout << " o!= 0. " << o << endl;
        return SIGN (o);
      }

    double a[3],b[3],c[3];
    int sign;
#if 0
    Point3 points[4];
    Vector3d p[4];
      
    points[0] = point1; 
    points[1] = point2; 
    points[2] = point3; 
    points[3] = point4;
    sign = sortPointIds(points, 4);

    p[0] = points[0].getVector();
    p[1] = points[1].getVector();
    p[2] = points[2].getVector();
    p[3] = points[3].getVector();
#else // point3 can directly access coordinates so no need for vector
    Point3 p[4];

    p[0] = point1; 
    p[1] = point2; 
    p[2] = point3; 
    p[3] = point4;
    sign = sortPointIds(p, 4);
#endif

    /* epsilon^1/8 */
    a[0] = p[1].position()(0); a[1] = p[1].position()(1);
    b[0] = p[2].position()(0); b[1] = p[2].position()(1);
    c[0] = p[3].position()(0); c[1] = p[3].position()(1);
    o = orient2d (a, b, c);
    //cout << "  /* epsilon^1/8 */ " << o << endl;
    if (o != 0.)
      return SIGN (o)*sign;


    /* epsilon^1/4 */
    a[0] = p[1].position()(0); a[1] = p[1].position()(2);
    b[0] = p[2].position()(0); b[1] = p[2].position()(2);
    c[0] = p[3].position()(0); c[1] = p[3].position()(2);
    o = orient2d (a, b, c);
    //cout << "   /* epsilon^1/4 */ " << o << endl;
    if (o != 0.)
      return - SIGN (o)*sign;
      
    /* epsilon^1/2 */
    a[0] = p[1].position()(1); a[1] = p[1].position()(2);
    b[0] = p[2].position()(1); b[1] = p[2].position()(2);
    c[0] = p[3].position()(1); c[1] = p[3].position()(2);
    o = orient2d (a, b, c);
    //cout << "   /* epsilon^1/2 */ " << o << endl;
    if (o != 0.)
      return SIGN (o)*sign;
      
    /* epsilon */
    a[0] = p[0].position()(0); a[1] = p[0].position()(1);
    b[0] = p[2].position()(0); b[1] = p[2].position()(1);
    c[0] = p[3].position()(0); c[1] = p[3].position()(1);
    o = orient2d (a, b, c);
    //cout << "   /* epsilon */ " << o << endl;
    if (o != 0.)
      return - SIGN (o)*sign;
      
    /* epsilon^5/4 */
    o = ORIENT1D (p[2].position()(0), p[3].position()(0));
    //cout << "    /* epsilon^5/4 */ " << o << endl;
    if (o != 0.)
      return SIGN (o)*sign;

    /* epsilon^3/2 */
    o = ORIENT1D (p[2].position()(1), p[3].position()(1));
    //cout << "    /* epsilon^3/2 */ " << o << endl;
    if (o != 0.)
      return - SIGN (o)*sign;
      
    /* epsilon^2 */
    a[0] = p[0].position()(0); a[1] = p[0].position()(2);
    b[0] = p[2].position()(0); b[1] = p[2].position()(2);
    c[0] = p[3].position()(0); c[1] = p[3].position()(2);
    o = orient2d (a, b, c);
    //cout << "    /* epsilon^2 */ " << o << endl;
    if (o != 0.)
      return SIGN (o)*sign;
      
    /* epsilon^5/2 */
    o = ORIENT1D (p[2].position()(2), p[3].position()(2));
    //cout << "     /* epsilon^5/2 */ " << o << endl;
    if (o != 0.)
      return SIGN (o)*sign;
      
    /* epsilon^4 */
    a[0] = p[0].position()(1); a[1] = p[0].position()(2);
    b[0] = p[2].position()(1); b[1] = p[2].position()(2);
    c[0] = p[3].position()(1); c[1] = p[3].position()(2);
    o = orient2d (a, b, c);
    //cout << "     /* epsilon^4 */ " << o << endl;
    if (o != 0.)
      return - SIGN (o)*sign;
      
    /* epsilon^8 */
    a[0] = p[0].position()(0); a[1] = p[0].position()(1);
    b[0] = p[1].position()(0); b[1] = p[1].position()(1);
    c[0] = p[3].position()(0); c[1] = p[3].position()(1);
    o = orient2d (a, b, c);
    //cout << "      /* epsilon^8 */ " << o << endl;
    if (o != 0.)
      return SIGN (o)*sign;
      
    /* epsilon^33/4 */
    o = ORIENT1D (p[1].position()(0), p[3].position()(0));
    //cout << "      /* epsilon^33/4 *  " << o << endl;
    if (o != 0.)
      return - SIGN (o)*sign;
      
    /* epsilon^17/2 */
    o = ORIENT1D (p[1].position()(1), p[3].position()(1));
    //cout << "     /* epsilon^17/2 */  " << o << endl;
    if (o != 0.)
      return SIGN (o)*sign;
      
    /* epsilon^10 */
    o = ORIENT1D (p[0].position()(0), p[3].position()(0));
    //cout << "      /* epsilon^10 *  " << o << endl;
    if (o != 0.)
      return SIGN (o)*sign;
      
    /* epsilon^21/2 */
    //cout << "     /* epsilon^21/2 */  " << sign <<  endl;
    return sign;
  }
}
