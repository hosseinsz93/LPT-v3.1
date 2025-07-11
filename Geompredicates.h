#ifndef _GEOM_PREDICATES_H
#define _GEOM_PREDICATES_H

namespace Geom
{
	double orient2d(double* pa, double* pb, double* pc);
	double orient3d(double* pa, double* pb, double* pc, double* pd);
	double incircle(double* pa, double* pb, double* pc, double* pd);
	double insphere(double* pa, double* pb, double* pc, double* pd, double* pe);

	void exactinit(void);
}

#endif


/* Automatic setting of emacs local variables. */
/* Local Variables: */
/* tab-width:2 */
/* End: */
