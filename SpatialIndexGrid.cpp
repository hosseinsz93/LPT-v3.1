#include "SpatialIndexGrid.h"

#include <limits>

#include "variables.h"

#include "Id.h"

#include "GeomCoordSystem.h"
#include "GeomFaceCalculations.h"
#include "GeomXFaceCalculations.h"
#include "GeomUtility.h"

#include "Ray3Track.h"
#include "RayFaceIntersect.h"

constexpr int SpatialIndexGrid::f[6][4];
constexpr int SpatialIndexGrid::n[8][3];
constexpr int SpatialIndexGrid::e[12][2][3];

int SpatialIndexGrid::ni = 0;
int SpatialIndexGrid::nj = 0;
int SpatialIndexGrid::nk = 0;
int SpatialIndexGrid::xs = 0;
int SpatialIndexGrid::xe = 0;
int SpatialIndexGrid::ys = 0;
int SpatialIndexGrid::ye = 0;
int SpatialIndexGrid::zs = 0;
int SpatialIndexGrid::ze = 0;
int SpatialIndexGrid::mx = 0;
int SpatialIndexGrid::my = 0;
int SpatialIndexGrid::mz = 0;

SpatialIndexGrid::SpatialIndexGrid(std::shared_ptr< UserCtx > & _user,
                                   std::shared_ptr< SpatialIndexGrid > & _self)
  : user(_user)
  , coor(user->coor)
{
  self.reset(this);
  _self = self;
  
  ni = user->IM;
  nj = user->JM;
  nk = user->KM;

  Id::setDim(ni, nj);

  DALocalInfo const & info = user->info;
  xs = info.xs, xe = info.xs + info.xm;
  ys = info.ys, ye = info.ys + info.ym;
  zs = info.zs, ze = info.zs + info.zm;
  mx = info.mx, my = info.my, mz = info.mz;

  maxDelta = std::max(std::abs(coor[0][0][0].x - coor[1][0][0].x),
                      std::abs(coor[0][0][0].y - coor[0][0][1].y));
  maxDelta = std::max(maxDelta, std::abs(coor[0][0][0].z - coor[0][1][0].z));
}

std::shared_ptr<Geom::SpatialIndexSource> SpatialIndexGrid::clone() const
{
  return self;
}

bool SpatialIndexGrid::isPoint(Geom::Id const &) const
{
  return false;
}

Vector3d SpatialIndexGrid::getCoordinate(Geom::Id const &id) const
{
  return Vector3d(std::numeric_limits<double>::max());
}

Geom::BoundingBox SpatialIndexGrid::getBoundingBox(Geom::Id const &id) const
{
  int i, j, k;
  Id(id).ijk(i, j, k);

  Cmpnts const & c = coor[k][j][i];
  Geom::BoundingBox bb(Vector3d(c.x, c.y, c.z));

  for (int idx=1; idx<8; ++idx)
  {
    Cmpnts const & c = coor [k+n[idx][2]] [j+n[idx][1]] [i+n[idx][0]];
    bb.expand(Vector3d(c.x, c.y, c.z));
  }

  return bb;
}

double SpatialIndexGrid::getDistanceSquared(Geom::Id const &id,
                                            Vector3d const &p) const
{
  // Doing this tricky assignments so that i,j,k can't accidentally
  // be used as loop indices.
  int ii, jj, kk;
  Id(id).ijk(ii, jj, kk);
  int const i = ii;
  int const j = jj;
  int const k = kk;

  int outsideCnt = 0;
  GeomX::Sidedness side[6];
  Vector3d area, center;
  std::vector<Vector3d> fcoor(4);
  for (int ii=0; ii<6; ++ii)
  {
    getFaceCoords(ii, i, j, k, user->coor, fcoor);
    faceAreaCentroid(fcoor, area, center, 1.0e-15);
    Geom::CoordSystem cs(center, area);
    auto pcs = cs.transInto(p);
    for (auto & f : fcoor)
      f = cs.transInto(f);

    if (pcs(2) > 1.0e-15)
    {
      side[ii] = GeomX::OUTSIDE;
      ++outsideCnt;
    }
    else
      side[ii] = GeomX::UNKNOWN;
  }
  if (outsideCnt == 0)
    return 0.0;

  // At this point, the point is outside the cell.  Time to find the distance.
  // Only need to find distances for faces where the point is marked as outside
  // the cell.
  double minDist2 = std::numeric_limits<double>::max();

  // We only need to look at faces that where tagged as outside since
  // those are the only ones the point sees.
  for (int ii=0; ii<6; ++ii)
  {
    if (side[ii] == GeomX::OUTSIDE)
    {
      // If the point is inside the face in a 2D sense using the coordinate
      // system defined by the face area normal and centroid, the distance
      // to the face is the closest.  Otherwise, the edges need to be checked
      // (which will also check the vertices).
      getFaceCoords(ii, i, j, k, user->coor, fcoor);
      faceAreaCentroid(fcoor, area, center, 1.0e-15);
      Geom::CoordSystem cs(center, area);
      auto pcs = cs.transInto(p);
      for (auto & f : fcoor)
        f = cs.transInto(f);

      double dist2 = std::numeric_limits<double>::max();
      auto fside = GeomX::pointSideofFace(fcoor, pcs, 1.0e-15);
      if (fside == GeomX::OUTSIDE)
      {
        for (int ee=0; ee<4; ++ee)
        {
          auto ldist2 = Geom::edgePointDistanceSquared(
                                              fcoor[ee], fcoor[(ee+1)%4], pcs);
          dist2 = std::min(dist2, ldist2);
        }
      }
      else
        dist2 = Geom::facePointDistanceSquared(fcoor, pcs);

      minDist2 = std::min(minDist2, dist2);
    }
  }

  return minDist2;
}

double SpatialIndexGrid::getDistanceSquaredId(Geom::Id const & id,
                                              Geom::Id queryId) const
{
  return std::numeric_limits<double>::max();
}

bool SpatialIndexGrid::intersectsRay(Geom::Id const & id,
                                     Geom::Ray3 const & _ray,
                                     double * argDist) const
{
  Ray3Track & ray = dynamic_cast<Ray3Track &>(const_cast<Geom::Ray3 &>(_ray));

  // Need the distance here so that we can get the face closest to the
  // ray's origin.  If there isn't one, make one.
  double locDist;
  double * dist = argDist ? argDist : &locDist;
  *dist = std::numeric_limits<double>::max();

  // Doing this tricky assignments so that i,j,k can't accidentally
  // be used as loop indices.
  int ii, jj, kk;
  Id(id).ijk(ii, jj, kk);
  int const i = ii;
  int const j = jj;
  int const k = kk;

  // Look at all the faces for this cell and see if the ray intersects
  // any of them.  If yes, choose the one that closest to the ray's origin.
  ray.setFaceId(-1);
  std::vector<Vector3d> fcoor(4);
  for (int ii=0; ii<6; ++ii)
  {
    getFaceCoords(ii, i, j, k, user->coor, fcoor);
    double wkDist;
    auto retVal = rayFaceIntersect(ray, fcoor, &wkDist);
    if (retVal && wkDist < *dist)
    {
      ray.setFaceId(ii);
      *dist = wkDist;
    }
  }

  return *dist < std::numeric_limits<double>::max();
}


void SpatialIndexGrid::getCellCoords(int const i, int const j, int const k,
                                     Cmpnts const * const * const * const coor,
                                     std::vector<Vector3d> & fcoor)
{
  constexpr int nVerts = sizeof(n) / sizeof(n[0]);
  fcoor.resize(nVerts);
  for (int ii=0; ii<nVerts; ++ii)
  {
    auto const & c = coor [k + n[ii][2]] [j + n[ii][1]] [i + n[ii][0]];
    fcoor[ii] = Vector3d(c.x, c.y, c.z);
  }
}

void SpatialIndexGrid::getFaceCoords(int faceId,
                                     int const i, int const j, int const k,
                                     Cmpnts const * const * const * const coor,
                                     std::vector<Vector3d> & fcoor)
{
    int const (&faceIdx)[4] = f[faceId];
    for (int ii=0; ii<4; ++ii)
    {
      // II used NI here to not interfere with class variable.
      int const (&NI)[3] = n[faceIdx[ii]];
      auto const & c = coor [k+NI[2]] [j+NI[1]] [i+NI[0]];
      fcoor[ii] = Vector3d(c.x, c.y, c.z);
    }
}

void SpatialIndexGrid::getVertCoords(int vertId,
                                     int const i, int const j, int const k,
                                     Cmpnts const * const * const * const coor,
                                     Vector3d & vcoor)
{
  // II used NI here to not interfere with class variable.
  int const (&NI)[3] = n[vertId];
  auto const & c = coor [k+NI[2]] [j+NI[1]] [i+NI[0]];
  vcoor = Vector3d(c.x, c.y, c.z);
}

SpatialIndexGrid::Sidedness SpatialIndexGrid::getCellSidedness(
                                              int i, int j, int k, int & iface)
{
  iface = -1;
  
  if (i == 0)
  {
    iface = 1;
    return Sidedness(ONBOUNDARY);
  }
  
  if (i == ni-2)
  {
    iface = 0;
    return Sidedness(ONBOUNDARY);
  }
  if (j == 0)
  {
    iface = 3;
    return Sidedness(ONBOUNDARY);
  }
  if (j == nj-2)
  {
    iface = 2;
    return Sidedness(ONBOUNDARY);
  }
  if (k == 0)
  {
    iface = 5;
    return Sidedness(ONBOUNDARY);
  }
  if (k == nk-2)
  {
    iface = 4;
    return Sidedness(ONBOUNDARY);
  }

  return Sidedness(INTERNAL);
}
