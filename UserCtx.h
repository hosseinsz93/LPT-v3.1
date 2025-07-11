#ifndef _USERCTX_H_
#define _USERCTX_H_

#include <limits>
#include <unordered_map>

/* This program is intended only for single CPU */

// From merge with dataPartExtr.cpp

struct BCS {
  Vec   Ubcs; // An ugly hack, waste of memory                                    
};

struct FlowWave
{
  PetscReal t, f;
};

class CrossX
{
  public:
    CrossX()
      : flag(0)
      , xmin(std::numeric_limits<double>::max())
      , ymin(std::numeric_limits<double>::max())
      , zmin(std::numeric_limits<double>::max())
      , xmax(std::numeric_limits<double>::min())
      , ymax(std::numeric_limits<double>::min())
      , zmax(std::numeric_limits<double>::min())
    { }
    
    void addMinLoc(double x, double y, double z);
    void addMaxLoc(double x, double y, double z);
    void makeInactive();

    void interpolate(double xTar, double & x, double & y, double & z);
    
    int const Flag() { return flag; }

    void const Min(double & x, double & y, double & z)
    {
      x = xmin;
      y = ymin;
      z = zmin;
    }

    void const Max(double & x, double & y, double & z)
    {
      x = xmax;
      y = ymax;
      z = zmax;
    }
  
  private:
    int flag;   // 0: nothing has happened to this particle so far.
                // 1: a minimum x value has been found.
                // 2: a maximum x value has been found
                //    - don't do anything else to this particle.
                // 3: particle is no longer active.
    double xmin, ymin, zmin;
    double xmax, ymax, zmax;
};

// end  merge with dataPartExtr.cpp

/* This program is intended only for single CPU */
struct UserCtx
{
  FILE * rkdebug;

  DA da, fda;
  DALocalInfo info;

  PetscReal dt;     // dt from the control file - time elapsed for each timestep.
  PetscReal tsdt;   //  time elapsed for each ts timesteps.

  PetscInt tis, tie, ts;

  Vec Ucat;         // Cartesian velocity at cell centers
  Vec Ucat_c;       // Cartesian velocity at grid nodes
  Vec Ucat_o;       // Previous timestep Cartesian velocity at cell centers

  Vec Cent;         // Coordinates of cell center
  Vec Rad;

  // This variable is only used in dataPartExtr.cpp.
  // Currently, it should not be used in the tracking program
  Vec Nvert;

  Cmpnts ***coor;
  Cmpnts *** ucat;
  Cmpnts *** ucat_o;
  Cmpnts ***cent;
  PetscScalar ***rad;

  PetscInt IM, JM, KM;

  PetscInt ptc_num;
  PetscInt kplane;

  List *CellList;

  // Do a restart at this time step if > 0.0.  WRO 2021-03-16
  PetscReal ptc_rstart;
  
  // New velocity calculation  WRO 2021-01-07
  // Current default is 1
  // 1: use minimum of Cd and kinatic gas models.
  // 3: kinetic gas
  // 4: use Runge Kutta 4 method
  // Anything else: massless method (particles exactly follow flow, no inertial effects)
  PetscInt calcMethod;

  std::shared_ptr<SpatialIndexGrid> spatialIndexGrid;
  std::unique_ptr<Geom::Rtree> rtreeGrid;

  std::shared_ptr<Surfaces> surfs;
  std::shared_ptr<SpatialIndexSurfs> spatialIndexSurfs;
  std::unique_ptr<Geom::Rtree> rtreeSurfs;

  Util::Timer globTimer;

  Vector3d DAccGrav; // Acceration of gravity

  // Fluid properties
  double DDenFd;   // Fluid density
  double DVisFd;   // Fluid dynamic viscosity

  // Default particle properties
  double DRoPt;    // particle density
  double DDiaPt;   // particle diameter

  // Scaling length and velocity for non-dimensional runs.
  double lengScale;
  double velScale;
    
  // Include ib surfaces as barriers for particle motion.
  double useIBs;

  // Method used to initialize particles to be tracked.
  // 0: Don't read ParticleInitial.dat
  // 1: Read from ParticleInitial.dat (default)
  double initPart;

  // Use methods to generate/modify particle infomation for Ali's covit
  // breath simulations.  Current, particles are not generated.
  // 1: Use breath routines
  double useBreath;

  // If this option is set, the particles' velocity is set to the fluid
  // velocity its the location plus the terminal velocity due to gravity.
  double setPartVel;
    
  // Debug drag calculation
  double useDebugDrag;
  std::unique_ptr<DragCalcDebug> debugDrag;

  // Debug particle trajectory (add info to trajectory outout files.)
  double useDebugTraj;

  // Result file output frequency in CFD timestep "units"
  int output_ts;

  // Run the Runge-Kutta integration routines.
  int RKtest;

  // Keep track of specific errors to give notification at end of program.
  unsigned globalError;

  // If set, output .vtk files instead of tecplot files.
  int vtk;

// From merge with dataPartExtr.cpp
  DA fda2;

  // If set, calculate the number of particles that went out the outlet.
  PetscInt end;

  PetscReal ***count;
    
  Vec   Csi, Eta, Zet, Aj;
  Vec   ICsi, IEta, IZet, IAj;
  Vec   JCsi, JEta, JZet, JAj;
  Vec   KCsi, KEta, KZet, KAj;
  Vec   Ucont;  // Contravariant velocity components
  Vec   P;
  Vec   Phi;
  Vec   GridSpace;
  BCS   Bcs;

  PetscReal     ren;    // Reynolds number
  PetscReal     st;     // Strouhal number

  PetscReal     r[101], tin[101], uinr[101][1001];

  Vec   lUcont, lUcat, lP, lPhi;
  Vec   lCsi, lEta, lZet, lAj;
  Vec   lICsi, lIEta, lIZet, lIAj;
  Vec   lJCsi, lJEta, lJZet, lJAj;
  Vec   lKCsi, lKEta, lKZet, lKAj;
  Vec   lGridSpace;
  Vec   lNvert, lNvert_o;
  Vec   lCent;

  Vec Ucat_sum;         // u, v, w
  Vec Ucat_cross_sum;           // uv, vw, wu
  Vec Ucat_square_sum;  // u^2, v^2, w^2

  PetscInt _this;

  FlowWave *inflow, *kinematics;
  PetscInt number_flowwave, number_kinematics;

  // Option to show count of points in cells that went through a cross plane
  // perpendicular to the x axis.
  bool xflag;
  PetscReal x;
  int partfile;

  // Option to show all the points that went thourgh a plane
  // perpendicular to the x axis, flow direction.
  bool xpartflag;
  PetscReal xpart;

  // Option to show all the points projected on a plane perpendicular to
  // the y axis for one time step for yfilesum or all files for ytitlsum.
  // (what I generally call the side view)
  bool yconc, yavrg;
  PetscReal y;

  // Option to show all the points projected on a plane perpendicular to
  // the z axis for one time step for zfilesum or all files for ztitlsum.
  // (what I generally call the top view)
  // For ground option:
  // For zavrg, only include particles that are close to the ground to get
  // the concentration on the ground.
  bool zconc, zavrg, groundFlag;
  PetscReal z, ground;

  // Scale for avrg options. NOTE: program also scales by dividing by the number
  // file. Scale is in addition to this.
  bool scaleFlag;
  PetscReal scale;

  // Option to show all the points with a z value < zmax. The points are
  // projected on a plane perpendicular to the z axis for one time step.
  // This option is used to find the concentrations a ground level.
  bool zmaxflag;
  PetscReal zmax;

  // Generate results file with only nvfield and ufield files.
  // As include FTLE results it ftle.output file which contains
  //   the ftle filenames to include in results.
  int uonly;

  std::unordered_map < int, CrossX > partMap;

// end  merge with dataPartExtr.cpp
};

#endif
