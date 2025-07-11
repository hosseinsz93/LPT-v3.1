#include "SpatialIndexBB.h"

#include <limits>

#include "GeomFaceCalculations.h"
#include "GeomXTol.h"

#include "Id.h"
#include "Surface.h"


template < typename T > T squ(T const & a) { return a*a; }

SpatialIndexBB::SpatialIndexBB(std::shared_ptr<Surfaces> & _surfs,
                               std::shared_ptr< SpatialIndexBB> & _self)
  : surfs(_surfs)
{
  self.reset(this);
  _self = self;
}

std::shared_ptr<Geom::SpatialIndexSource> SpatialIndexBB::clone() const
{
  return self;
}

bool SpatialIndexBB::isPoint(Geom::Id const &) const
{
  return false;
}

Vector3d SpatialIndexBB::getCoordinate(Geom::Id const &id) const
{
  return Vector3d(std::numeric_limits<double>::max());
}

Geom::BoundingBox SpatialIndexBB::getBoundingBox(Geom::Id const &id) const
{
  int zId = id.getValue() - 1;
  return surfs->getSurf(zId)->getbb();
}

double SpatialIndexBB::getDistanceSquared(Geom::Id const &id,
                                          Vector3d const &p) const
{
  int zId = id.getValue() - 1;
  auto bb = surfs->getSurf(zId)->getbb();

  double dist2 = 0.0;

  if (p(0) < GeomX::Tol(bb.min()(0)))
    dist2 += squ(p(0) - bb.min()(0));
  if (p(0) > GeomX::Tol(bb.max()(0)))
    dist2 += squ(p(0) - bb.max()(0));

  if (p(1) < GeomX::Tol(bb.min()(1)))
    dist2 += squ(p(1) - bb.min()(1));
  if (p(1) > GeomX::Tol(bb.max()(1)))
    dist2 += squ(p(1) - bb.max()(1));

  if (p(2) < GeomX::Tol(bb.min()(2)))
    dist2 += squ(p(2) - bb.min()(2));
  if (p(2) > GeomX::Tol(bb.max()(2)))
    dist2 += squ(p(2) - bb.max()(2));

  return dist2;
}

bool SpatialIndexBB::intersectsRay(Geom::Id const &/*id*/,
                                   Geom::Ray3 const & /*ray*/,
                                   double * /*distance*/) const
{
  return false;
}
