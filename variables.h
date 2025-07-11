#ifndef _VARIABLES_H
#define _VARIABLES_H

#include "petscvec.h"
#include "petscda.h"
#include "petscksp.h"
#include "petscsnes.h"

#include <stdio.h>

#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "timer.h"
#include "cmpnts.h"
#include "GeomVector.h"

// Current version of program.
constexpr char const * progVersion = "20221221a";

// Predefine these classes to avoid circular includes.
namespace Geom
{
  class Rtree;
}
class Surfaces;
class SpatialIndexGrid;
class SpatialIndexSurfs;
class DragCalcDebug;

//#include "petscvec.h"
//#include "petscda.h"
//#include "petscksp.h"
//#include "petscsnes.h"

enum GeomLookupStatus
{
  GEOM_NOT_DONE,
  GEOM_SUCCESSFUL, // The closes cell has been found and the i,j,k is set.
  GEOM_FAILED      // Geom Lookup failed to find a cell close enough.
};

enum ParticleStatus
{
  INSIDE_ACTIVE,         //  0
  LOOKUP_FAILED,         //  1
  ZERO_VEL_LOC,          //  2
  IN_IB_CELL,            //  3
  IN_SOLID_CELL,         //  4
  OUTSIDE_AT_INLET,      //  5
  OUTSIDE_AT_WALL,       //  6
  OUTSIDE_AT_OUTLET,     //  7
  ON_IB_SURFACE,         //  8
  OUTSIDE_UNKNOWN,       //  9
  ZERO_DIAMETER          // 10
};

static constexpr char const *partMess[]
{
  "inside active",
  "lookup failed",
  "zero velocity",
  "in an ib cell",
  "in a solidcell",
  "outside at inlet",
  "outside at wall",
  "outside at outlet",
  "on ib surface",
  "outside unknown",
  "zero diameter"
};

enum GlobalErrors
{
  MOVED_TOO_FAR         = 0x01,         // At least one particle was moved too
                                        // far. Consider reducing the time
                                        // step.
  GRID_QUERY_RAYTRACE2  = 0x02,         // Grid query in raytrace2 failed.
  RAY_QUERY_RAYTRACE2   = 0x04,         // Ray query in raytrace2 failed.
  FACE_ID_RAYTRACE2     = 0x08,         // Problem retreiving face Id in
                                        // raytrace2.
  SURF_INTERS_PART_TRAJ = 0x10,         // Problem with surface intersection in
                                        // ParticleTrajectory.
  IB_GRID_LOOKUP_ERROR  = 0x20,         // Error in grid lookup in ib routine.
  RAY_QUERY_PART_TRAJ   = 0x40,         // Ray query error in
                                        // ParticleTrajectory.
  COEF_ERROR            = 0x80,         // Error of coef calculation in
                                        // ISPointInCell
  FORCE_REL_VEC_TO_ZERO = 0x0100        // Because of a very small reynold's
                                        // number, the particle's relative
                                        // to the air flow was set to zero.
};

extern PetscInt binary_input;
extern PetscInt block_number;

// This is a base class used for polymorphism
class Object
{
};

class Vector3ds
{
  private:
    Vector3d x;

  public:
    Vector3ds() : x(0.0) { }
    Vector3ds(Vector3d  const & _x) : x(_x)   { }
    Vector3ds(Vector3ds const & _x) : x(_x.x) { }
    
    Vector3ds const & get(bool const debug = false) const
    {
      if (debug)
        int bug = 1;
      return *this;
    }
    
    Vector3d const & getv(bool const debug = false) const
    {
      if (debug)
        int bug = 1;
      return *this;
    }

    double operator()(int idx)
    {
      return x(idx);
    }

    operator Vector3d       & ()       { return x; }
    operator Vector3d const & () const { return x; }

    void set(Vector3d const & _x, bool const debug = false)
    {
      if (debug)
        int bug = 1;
      x = _x;
    }

    void set(Cmpnts const & _x, bool const debug = false)
    {
      if (debug)
        int bug = 1;
      x(0) = _x.x;
      x(1) = _x.y;
      x(2) = _x.z;
    }

    void set(double const _x, double const _y, double const _z,
             bool const debug = false)
    {
      if (debug)
        int bug = 1;
      x(0) = _x;
      x(1) = _y;
      x(2) = _z;
    }

    void set(int const idx, double const _x, bool const debug = false)
    {
      if (debug)
        int bug = 1;
      x(idx) = _x;
    }
};



struct Location
{
  Location()                       : i(-1), j(-1), k(-1) {}
  Location(int _i, int _j, int _k) : i(_i), j(_j), k(_k) {}
  PetscInt i, j, k;
};

struct node
{
  Location Node;
  struct node *next;
};

struct List
{
  node *head;
  PetscInt number;
};


struct Node_List
{
  PetscInt	index;
  struct list_node *next;
};


class Particle
{
  public:
  Particle(ParticleStatus _partStat = INSIDE_ACTIVE)
  {
    coor.x = coor.y = coor.z = 0.0;
    coef.x = coef.y = coef.z = 0.0;
    Rey = 0.0;
    RespTime = 0.0;
    dt = 0.0;
    next_time = 0;
    firstVel = true;
    stime = 0;
    etime = 0;
    DDiaPt = 0;
    DDenPt = 0;
    dist = 0;
    partStat = _partStat;
    geomLookupStatus = GEOM_NOT_DONE;
  }

  Cmpnts coor;
  Location loc;
  Cmpnts coef;
  Vector3ds Vf0;   // Fluid velocity - previous time step
  Vector3ds Vf1;   // Fluid velocity - current time step
  Vector3ds Vp;    // Particle velocity
  double Rey;      // Reynold's number at beginning of iteration.
  double RespTime; // Particle's response time at beginning of iteration.

  // I'm turning this off for now to see if the program doesn't blow up anymore.
  // std::vector < Cmpnts > ztraj;
  PetscReal dt; /* time step used to integrate the particle trajectory */
  PetscInt next_time;
  PetscInt firstVel;

  // Particle property variables I added.  WRO 2021-01-06
  PetscReal stime;
  PetscReal etime;
  PetscReal DDiaPt;
  PetscReal DDenPt;

  // The distance returned by an successful call to rtreeGrid->getClosestEntity
  PetscReal dist;
    
  // INSIDE_ACTIVE,
  // ZERO_VEL_LOC,        point is at a zero velocity location
  // IN_IB_CELL,          point is in an IB cell
  // IN_SOLID_CELL,       point is in a solid cell
  // OUTSIDE_AT_WALL,     point is outside the domain - wall
  // OUTSIDE_AT_OUTLET    point is outside the domain - outlet
  // ON_IB_SURFACE        point is on ib surface
  ParticleStatus partStat;

  // These values must be set correctly before going into ParticleLocation
  // subroutine.  The show the status of the geometry lookup to which finds
  // the i,j,k value for the particles current location.
  // Allowable values are defined in enum GeomLookupStatus.
  // GEOM_NOT_DONE,
  // GEOM_SUCCESSFUL, The closes cell has been found and the i,j,k is set.
  // GEOM_FAILED      Geom Lookup failed to find a cell close enough.
  GeomLookupStatus geomLookupStatus;


  std::string toStr()
  {
    std::stringstream sout;
    sout << "coor:  " << coor.x << "," << coor.y << "," << coor.z << "\n";
    sout << "loc:   " << loc.i  << "," << loc.j  << "," << loc.k  << "\n";
    sout << "coef:  " << coef.x << "," << coef.y << "," << coef.z << "\n";
    sout << "Vf0:   " << Vf0.getv()(0) << "," << Vf0.getv()(1) << "," << Vf0.getv()(2) << "\n";
    sout << "Vf1    " << Vf1.getv()(0) << "," << Vf1.getv()(1) << "," << Vf1.getv()(2) << "\n";
    sout << "Vp:    " << Vp.getv()(0) << "," << Vp.getv()(1) << "," << Vp.getv()(2) << "\n";
    sout << "dt:    " << dt     << "\n";
    sout << "next_time: " << next_time << "\n";
    sout << "stime: " << stime << "\n";
    sout << "etime: " << etime << "\n";
    sout << "DDiaPt: " << DDiaPt << "\n";
    sout << "DDenPt: " << DDenPt << "\n";
    sout << "partStat: " << partStat << ": " << partMess[partStat] << "\n";
    return sout.str();
  }
};

#include "UserCtx.h"

template < typename T >
T CENT(T const * const * const * const var, int k, int j, int i)
{
  return ((var[k  ][j  ][i  ] +
           var[k  ][j  ][i+1] +
           var[k  ][j+1][i  ] +
           var[k  ][j+1][i+1] +
           var[k+1][j  ][i  ] +
           var[k+1][j  ][i+1] +
           var[k+1][j+1][i  ] +
           var[k+1][j+1][i+1]) * 0.125);
}


// Deleter for std::shared_ptr.  This is needed since petsc allocates
// some memory that is put into a shared pointer we don't want the shared
// pointer deleting it.  In other words, the shared pointer is not the
// owner of the memory.
template < typename T > struct NoDeleter
{ 
  void operator()(T* p) const { /* Do nothing! */ }
};

// Subroutine definitions.
PetscErrorCode Initial(FILE *fd, std::shared_ptr< UserCtx > & user,
                       Util::Timer & timer);
PetscErrorCode ParticleInitial(std::shared_ptr< UserCtx > & user,
                               std::vector<Particle> & Ptc);
PetscErrorCode ParticleLocation1(std::shared_ptr< UserCtx > & user,
                                 Particle *P, int pId);
PetscErrorCode VelocityReadin(std::shared_ptr< UserCtx > & user, PetscInt ti,
                              PetscInt bi);
PetscErrorCode tracking(std::shared_ptr< UserCtx > & user,
                        std::vector<Particle> & Ptc);


void initlist(List *ilist);
void insertnode(List *ilist, Location Node);

#endif
