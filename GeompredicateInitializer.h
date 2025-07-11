#ifndef _GEOM_PREDICATEINITIALIZER_H
#define _GEOM_PREDICATEINITIALIZER_H

#include "Geompredicates.h"

/* Need to call exactinit once per program instantiation for the robust
 * intersection stuff to work correctly.  Look at comments in predicates.cpp
 * for details.
 */

namespace Geom
{
  class PredicateInitializer
    {
    public:
      PredicateInitializer() {exactinit();}
    };

  static PredicateInitializer _initialize;
}

#endif


/* Automatic setting of emacs local variables. */
/* Local Variables: */
/* tab-width:2 */
/* End: */
