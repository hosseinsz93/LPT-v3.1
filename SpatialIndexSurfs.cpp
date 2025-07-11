#include "SpatialIndexSurfs.h"

#include <limits>

#include "Id.h"

#include "variables.h"
#include "GeomCoordSystem.h"
#include "GeomFaceCalculations.h"
#include "GeomXFaceCalculations.h"
#include "GeomUtility.h"
#include "RayFaceIntersect.h"
#include "SpatialIndexGrid.h"


SpatialIndexSurfs::SpatialIndexSurfs(std::shared_ptr<UserCtx> & _user,
                                     std::shared_ptr<SpatialIndexSurfs> & _self)
  : user(_user)
  , surfs(_user->surfs)
{
  self.reset(this);
  _self = self;
}

std::shared_ptr<Geom::SpatialIndexSource> SpatialIndexSurfs::clone() const
{
  return self;
}

bool SpatialIndexSurfs::isPoint(Geom::Id const &) const
{
  return false;
}

Vector3d SpatialIndexSurfs::getCoordinate(Geom::Id const &id) const
{
  return Vector3d(std::numeric_limits<double>::max());
}

Geom::BoundingBox SpatialIndexSurfs::getBoundingBox(Geom::Id const &id) const
{
  int sId=0, zId=0;
  Geom::BoundingBox bb;
  std::tie(sId, zId) = surfs->convTId(id);
  auto coors = surfs->getSurf(sId)->getTriCoor(zId);
  bb.expand(coors[0]);
  bb.expand(coors[1]);
  bb.expand(coors[2]);
  return bb;
}

double SpatialIndexSurfs::getDistanceSquared(Geom::Id const &id,
                                             Vector3d const &p) const
{
  int sId=0, zId=0;
  std::tie(sId, zId) = surfs->convTId(id);
  auto coors = surfs->getSurf(sId)->getTriCoor(zId);
  auto dist2 = Geom::trianglePointDistanceSquared(coors[0],coors[1],coors[2],p);

#if 0
  PetscPrintf(PETSC_COMM_WORLD, "\nid: %d  sId: %d  zId: %d  dist: %f\n",
              id.getValue(), sId, zId, sqrt(dist2));
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

// Don't define this function since queryId is not defined or used
// in the system of this program.
// double SpatialIndexSurfs::getDistanceSquaredId(Geom::Id const & id,
//                                                Geom::Id queryId) const;


bool SpatialIndexSurfs::intersectsRay(Geom::Id const & id,
                                      Geom::Ray3 const & ray,
                                      double * dist) const
{
  int sId=0, zId=0;
  std::tie(sId, zId) = surfs->convTId(id);
  auto coors = surfs->getSurf(sId)->getTriCoor(zId);
  auto retVal = rayFaceIntersect(ray, coors, dist);
  if (!retVal)
  {
    *dist = std::numeric_limits<double>::max();
    return retVal;
  }

#if 0
  double searchRad = 3.0 * user->spatialIndexGrid->getMaxDelta() / 5.0;
  Vector3d area, centroid;
  Geom::triangleAreaCentroid(coors[0], coors[1], coors[2], area, centroid);
  area.normalize();

  PetscPrintf(PETSC_COMM_WORLD, "\nSpatialIndexSurfs::intersectsRay\n");
  PetscPrintf(PETSC_COMM_WORLD, "dist: %f\n", *dist);
  PetscPrintf(PETSC_COMM_WORLD, "view,%f,%f,%f\n", area(0), area(1), area(2));

  PetscPrintf(PETSC_COMM_WORLD, "v,1,%f,%f,%f\n",
              coors[0](0), coors[0](1), coors[0](2));
  PetscPrintf(PETSC_COMM_WORLD, "v,2,%f,%f,%f\n",
              coors[1](0), coors[1](1), coors[1](2));
  PetscPrintf(PETSC_COMM_WORLD, "v,3,%f,%f,%f\n",
              coors[2](0), coors[2](1), coors[2](2));

  PetscPrintf(PETSC_COMM_WORLD, "v,4,%f,%f,%f\n",
                                        centroid(0), centroid(1), centroid(2));
  auto tmp = centroid + area * searchRad;
  PetscPrintf(PETSC_COMM_WORLD, "v,5,%f,%f,%f\n", tmp(0), tmp(1), tmp(2));

  auto orig = ray.getOrigin();
  PetscPrintf(PETSC_COMM_WORLD, "v,6,%f,%f,%f\n", orig(0), orig(1), orig(2));
  tmp = orig + ray.getDirection() * searchRad;
  PetscPrintf(PETSC_COMM_WORLD, "v,7,%f,%f,%f\n", tmp(0), tmp(1), tmp(2));

  // Calculate the new location.
  tmp = orig + ray.getDirection() * *dist;
  PetscPrintf(PETSC_COMM_WORLD, "v,8,%f,%f,%f\n", tmp(0), tmp(1), tmp(2));
#endif

  return retVal;
}
