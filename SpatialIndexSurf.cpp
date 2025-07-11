#include "SpatialIndexSurf.h"

#include <limits>

#include "Id.h"

#include "GeomFaceCalculations.h"
#include "RayFaceIntersect.h"
#include "Surface.h"


SpatialIndexSurf::SpatialIndexSurf(std::shared_ptr<Surface> & _surf,
                                   std::shared_ptr< SpatialIndexSurf > & _self)
  : surf(_surf)
{
  self.reset(this);
  _self = self;
}

std::shared_ptr<Geom::SpatialIndexSource> SpatialIndexSurf::clone() const
{
  return self;
}

bool SpatialIndexSurf::isPoint(Geom::Id const &) const
{
  return false;
}

Vector3d SpatialIndexSurf::getCoordinate(Geom::Id const &id) const
{
  return Vector3d(std::numeric_limits<double>::max());
}

Geom::BoundingBox SpatialIndexSurf::getBoundingBox(Geom::Id const &id) const
{
  Geom::BoundingBox bb;
  auto const & t = surf->gettri(id);
  bb.expand(surf->getvert(t[0]));
  bb.expand(surf->getvert(t[1]));
  bb.expand(surf->getvert(t[2]));
  return bb;
}

double SpatialIndexSurf::getDistanceSquared(Geom::Id const &id,
                                            Vector3d const &p) const
{
  auto coors = surf->getTriCoor(id);
  auto dist2 = Geom::trianglePointDistanceSquared(coors[0],coors[1],coors[2],p);

#if 0
  PetscPrintf(PETSC_COMM_WORLD, "\ndist: %f\n", sqrt(dist2));
  PetscPrintf(PETSC_COMM_WORLD, "v,1,%f,%f,%f\n", p(0), p(1), p(2));
  PetscPrintf(PETSC_COMM_WORLD, "v,2,%f,%f,%f\n",
              coors[0](0), coors[0](1), coors[0](2));
  PetscPrintf(PETSC_COMM_WORLD, "v,3,%f,%f,%f\n",
              coors[1](0), coors[1](1), coors[1](2));
  PetscPrintf(PETSC_COMM_WORLD, "v,4,%f,%f,%f\n",
              coors[2](0), coors[2](1), coors[2](2));
#endif

  return dist2;
}

bool SpatialIndexSurf::intersectsRay(Geom::Id const & id,
                                     Geom::Ray3 const & ray,
                                     double * dist) const
{
  auto coors = surf->getTriCoor(id);
  return rayFaceIntersect(ray, coors, dist);
}
