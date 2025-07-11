#ifndef _RAYFACEINTERSECT_H_
#define _RAYFACEINTERSECT_H_

#include <vector>

#include "GeomRay3.h"
#include "GeomVector.h"

bool rayFaceIntersect(Geom::Ray3 const & ray,
                      std::vector < Vector3d > & faceCoor,
                      double * dist = 0, Vector3d * intersect = 0);

#endif
