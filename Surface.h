#ifndef _SURFACE_H_
#define _SURFACE_H_

#include <array>
#include <fstream>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "GeomId.h"
#include "GeomBoundingBox.h"
#include "GeomVector.h"

struct UserCtx;
namespace Geom
{
  class Rtree;
}
class SpatialIndexBB;
class SpatialIndexSurf;

enum CoorLabel
{
  X, Y, Z
};


class Surface
{
  public:
    Surface(std::string _filename);

  private:
    void fileError();
    void readCoor(std::fstream & fin, CoorLabel coorIdx);
    void readNorm(std::fstream & fin);
    void readTris(std::fstream & fin);

  public:
    void createRtree(std::shared_ptr<Surface> & surf);

    unsigned vsize() const { return verts.size(); }
    unsigned tsize() const { return tris.size(); }

    std::vector < Vector3d > const & getverts() const { return verts; }
    Vector3d const & getvert(int vId) const { return verts[vId]; }

    std::vector < std::array<int,3> > const & gettris() const { return tris; }
    std::array<int,3> const & gettri(int tId) const { return tris[tId]; }
    std::array<int,3> const & gettri(Geom::Id tId) const
      { return gettri(tId.getValue() - 1); }
    
    static void remCargRtn(std::string & buff);
    void outputStar() const;

    Geom::BoundingBox const & getbb() const { return bb; }

    // Returns the coordinates of a triangle.
    std::vector < Vector3d > getTriCoor(Geom::Id const & id) const
                                      { return getTriCoor(id.getValue() - 1); }
    std::vector < Vector3d > getTriCoor(int tId) const;

    Geom::Rtree * getRtreeSurf() { return rtreeSurf.get(); }
    

  private:
    std::string filename;
    int nverts, ntris;
    std::vector < Vector3d > verts;
    std::vector < std::array<int,3> > tris;

    Geom::BoundingBox bb;

    std::shared_ptr<SpatialIndexSurf> spatialIndexSurf;
    std::unique_ptr<Geom::Rtree> rtreeSurf;
};

class Surfaces
{
  public:
    Surfaces(std::shared_ptr< UserCtx > & _user,
             std::shared_ptr< Surfaces > & _self);
    
    int getSurfsSize() const { return surfs.size(); }
    std::shared_ptr < Surface > const & getSurf(int idx) const { return surfs[idx]; }

    std::vector < int > const & getTIdxs() const { return tIdxs; }

    std::tuple<int,int> convTId(Geom::Id const &id) const;

    // Returns the coordinates of a triangle.
    std::vector < Vector3d > getTriCoor(Geom::Id const &id) const;
  
  private:
    void createGlobalRtree();
    void createSurfBBRtree();

  private:
    std::shared_ptr< UserCtx > user;
    std::shared_ptr< Surfaces > self;
    
    std::vector < std::shared_ptr < Surface > > surfs;

    std::vector < int > vIdxs;
    std::vector < int > tIdxs;

    std::shared_ptr<SpatialIndexBB> spatialIndexBB;
    std::unique_ptr<Geom::Rtree> rtreeBB;
};

#endif
