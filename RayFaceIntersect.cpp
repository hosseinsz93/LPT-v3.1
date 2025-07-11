#include "RayFaceIntersect.h"

#include "GeomCoordSystem.h"
#include "GeomXFaceCalculations.h"
#include "GeomXTol.h"

#include "Ray3Track.h"


// See if ray and face intersect.  If not, return false. If yes, return true
// and the distance to and location of the intersection point.

// THIS ROUTINE DOESN'T WORK IF ITS FACE OVERLAPS ITSELF!

// This is a helper routine.
bool rayFaceTriIntersect(Vector3d const & p,
                         std::vector < Vector3d > & faceCoor,
                         double * dist, Vector3d * intersect);

bool rayFaceIntersect(Geom::Ray3 const & _ray,
                      std::vector < Vector3d > & faceCoor,
                      double * dist, Vector3d * intersect)
{
  Ray3Track const & ray =  dynamic_cast<Ray3Track const &>(_ray);

  // Transform coordinates in from global to cs.
  std::vector < Vector3d > csFace(faceCoor.size());
  for (unsigned i=0; i<faceCoor.size(); ++i)
    csFace[i] = ray.getCS().transInto(faceCoor[i]);


  if (faceCoor.size() == 3)
  {
    auto status = rayFaceTriIntersect(Vector3d(0.0), csFace, dist, intersect);
    if (intersect)
      *intersect = ray.getCS().transOutOf(*intersect);
    return status;
  }
  else
  {
    std::vector < Vector3d > tri(3);
    tri[2] = Geom::faceCentroid(csFace);
    for (unsigned i=0; i<csFace.size(); ++i)
    {
      tri[0] = csFace[i];
      tri[1] = csFace[(i+1)%csFace.size()];
      bool ret = rayFaceTriIntersect(Vector3d(0.0), tri, dist, intersect);
      if (intersect)
        *intersect = ray.getCS().transOutOf(*intersect);
      if (ret)
        return ret;
    }
  }

  return false;
}


bool rayFaceTriIntersect(Vector3d const & p,
                         std::vector < Vector3d > & faceCoor,
                         double * dist, Vector3d * intersect)
{
  auto areaCoor = GeomX::triangleAreaCoords2D(p, &faceCoor[0]);

  for (int i=0; i<3; ++i)
    if (areaCoor(i) < GeomX::Tol(0.0))
      return false;

  double csz = areaCoor(0) * faceCoor[0](2)
             + areaCoor(1) * faceCoor[1](2)
             + areaCoor(2) * faceCoor[2](2);

  // Make sure the point is in the path of the ray and not "behind" it.
  if (csz < GeomX::Tol(0.0))
    return false;

  if (dist || intersect)
  {
    Vector3d interPoint(0.0, 0.0, csz);
    if (dist)
      *dist = (interPoint - p).mag();
    if (intersect)
      *intersect = interPoint;
  }
  
  return true;
}
