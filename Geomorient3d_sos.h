#ifndef _GEOM_ORIENT3D_SOS_H
#define _GEOM_ORIENT3D_SOS_H


namespace Geom
{
  class Point3;

  double orient3d_sos(Point3 const &p1,
                      Point3 const &p2,
                      Point3 const &p3,
                      Point3 const &p4,
                      double *orient3dresult = 0);
}

#endif


/* Automatic setting of emacs local variables. */
/* Local Variables: */
/* tab-width:2 */
/* End: */
