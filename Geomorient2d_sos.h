#ifndef _GEOM_ORIENT2D_SOS_H
#define _GEOM_ORIENT2D_SOS_H




namespace Geom
{
  class Point2;

  double orient2d_sos(Point2 const &p1,
                      Point2 const &p2,
                      Point2 const &p3,
                      double *orient2dresult = 0);
}

#endif


/* Automatic setting of emacs local variables. */
/* Local Variables: */
/* tab-width:2 */
/* End: */
