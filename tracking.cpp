#include "variables.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <memory>

#include "GeomFaceCalculations.h"
#include "GeomId.h"
#include "GeomRtree.h"
#include "GeomUtility.h"
#include "GeomVector.h"
#include "GeomXTol.h"
#include "GeomXUtilities.h"

#include "BreathParticles.h"
#include "calcParticleAcceleration.h"
#include "DragCalcDebug.h"
#include "Id.h"
#include "Ray3Track.h"
#include "ReysToCd.h"
#include "RungeKutta.h"
#include "SpatialIndexGrid.h"
#include "Surface.h"

// Uncomment this to get output from ParticleTrajectory subroutine.
// #define WRO_DEBUG
// #define WRO_DEBUG_PID 1

#if defined(WRO_DEBUG_PID)
  constexpr int wroDebugPID = WRO_DEBUG_PID;
#elif defined(WRO_DEBUG)
  constexpr int wroDebugPID = -1;
#else
  constexpr int wroDebugPID = -2;
#endif


// Subroutine definitions.
PetscErrorCode BoundingSphere(std::shared_ptr< UserCtx > & user);

PetscErrorCode UniformGrid(std::shared_ptr< UserCtx > & user);

PetscErrorCode ParticleTrajectoryIO2(std::shared_ptr< UserCtx > & user,
                                     std::vector<Particle> & Ptc,
                                     PetscInt timestep, double time);

// Write particle trajectory tracks in new format
// timestep is the index to be used for the filename.
// time is the time of ti.
PetscErrorCode ParticleTrajectoryIO3(std::shared_ptr< UserCtx > & user,
                                     std::vector<Particle> & Ptc,
                                     PetscInt timestep, double time);

PetscTruth ISPointInCell(Cmpnts p, UserCtx * user,
                         Location loc, Cmpnts *coef, Particle *P, int pId,
                         bool printErr = true);

PetscErrorCode ParticleLocation(UserCtx * user, Particle *P,
                                PetscInt iradial, int pId, bool outsideIsAnError = true);

PetscErrorCode VelocityInterpolation(UserCtx * user,
                                     Particle *P, PetscReal time,
                                     PetscReal time_o, PetscReal time_n, bool VF0update = true);

PetscErrorCode ParticleTrajectory(std::shared_ptr< UserCtx > & user,
                                  Particle * const P, PetscInt ti, int pId,
                                  double time, double time_o, double time_n);

PetscErrorCode PFaceDis(Cmpnts p, Cmpnts p0, Cmpnts p1, Cmpnts p2, Cmpnts p3,
                        PetscReal *dis);

bool raytrace(Vector3d const & P0, Vector3d & P1, double & dist,
              std::unique_ptr<Geom::Rtree> & rtreeSurfs, Geom::Id & result,
              std::shared_ptr< UserCtx > & user, int pId /* for debugging only */ );

void raytrace2(Vector3d const & P0, Vector3d & P1,
               int pId, Particle * P, std::shared_ptr< UserCtx > & user,
               Geom::Id & result, Location & loc, double searchRad);

void PrintDebugInfo01(bool header, Particle *P, int pId, int ti,
                      Vector3d const & DAcc, Vector3d const & Vf,
                      Vector3d const & Vpf1, Vector3d const & Vp1,
                      Vector3d const & P1);

void PrintDebugInfo02(int pId, Vector3d const & P0, double ptc_dt,
                      double DVolPt, double DMasPt, Vector3d const & DAccGrav,
                      Vector3d const & Vf, Vector3d const & Vp0,
                      Vector3d const & Vpf0, double rey, double Cd,
                      Vector3d const & DFdrag, Vector3d const & DFtot,
                      Vector3d const & DAcc, Vector3d const & DVel,
                      Vector3d const & Vpf1, double const & Vt,
                      Vector3d const & Vp1, Vector3d const & Pdel,
                      Vector3d const & P1);

void updatePartStatus(int pId, ParticleStatus partStat, Location const & loc,
                      Vector3d const & coor, Particle *P);


template < typename T >
double Length(T v)
{
  return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}


// Calculate the terminal velocity and see how it compares to Vpf1.
// Adjust as necessary.
double calcTermVel(Vector3d const & DAcc, double const U,
                   std::shared_ptr< UserCtx > & user, Particle * const P)
{
  return DAcc.mag() * GeomX::pow(P->DDiaPt, 2) * (P->DDenPt - user->DDenFd)
     / (U * 18.0 * user->DVisFd);
}

// Add the terminal velocity to Vpf.
Vector3d addTermVel(double const Vt, Vector3d const & DAcc,
                    Vector3d const & Vpf)
{
  // If there is no accelleration, the relative velocity must be zero.
  if (Vt == GeomX::Tol(0.0))
    return Vector3d(0.0);

  return Vpf + Vt * DAcc / DAcc.mag();
}


// Calculate drag using Kinetic Gas Theory
// Vpf0 MUST be the real velocity not the scaled velocity!!!
Vector3d dragKineticGas
(double DVisFd, double DDiaPt, Vector3d const & Vpf0)
{
  // Mean free path length for air at STP
  constexpr double lamda = 0.066e-6;

  // Cunningham slip correction parameter
  constexpr double A = 1.234;
  constexpr double B = 0.414;
  constexpr double E = 0.870;
  constexpr double AB = A + B;

  // Calculate base drag force then divide by the correction factor.
  // Formula taken from Wang, Drag force, diffustion coefficient and
  // mobility of small particles..., physical review 68, 2003, equ. 5
  // Constants taken from: B. E. Dahneke, J. Aerosol Sci. 4, 139 1972
  Vector3d DFdrag = -3.0 * Geom::PI * DVisFd * DDiaPt * Vpf0;
  double DKn = lamda / (DDiaPt / 2.0);
  DFdrag /= 1.0 + DKn * (A + B * exp(-E/DKn));
  return DFdrag;
};


void genRtree(std::shared_ptr< UserCtx > & user);

extern PetscReal ptc_dt;
extern PetscInt ptc_tn;
// Number of particle iterations between the iterations of the results files.
extern PetscInt ptc_ftn;
// Number of particle interations in each particle iteration.
extern PetscInt ptc_int;

//extern PetscReal cl;	/* mli */
  

PetscErrorCode tracking(std::shared_ptr< UserCtx > & user,
                        std::vector<Particle> & Ptc)
{
  PetscInt bi, blocknumber=1;
  PetscInt ti, tistart=0.0;
  PetscInt  itime=0,  itime_o=0,  itime_n=0;
  PetscReal time=0.0, time_o=0.0, time_n=0.0;
  PetscInt pn;

  Vec Coor;
  DAGetGhostedCoordinates(user->da, &Coor);
  DAVecGetArray(user->fda, Coor, &user->coor);

  DAVecGetArray(user->da,  user->Rad,  &user->rad);
  DAVecGetArray(user->fda, user->Cent, &user->cent);

  PetscPrintf(PETSC_COMM_WORLD, "blocknumber %d \n",blocknumber);

#if 0
  for (bi=0; bi < blocknumber; bi++)
  {
    BoundingSphere(&user[bi]);
    UniformGrid(&user[bi]);
  }
#endif
  
  DAVecGetArray(user->fda, user->Ucat,   &user->ucat);
  DAVecGetArray(user->fda, user->Ucat_o, &user->ucat_o);

  
  genRtree(user);

  // Setup information for particle location here since that information
  // is not added during particle read.
  for (pn=0; pn<user->ptc_num; pn++) 
    ParticleLocation(user.get(), &Ptc[pn], 5, pn);

  // This section describes all the variables involved with time and
  // time stepping.  It is important to differentiate between time steps
  // of the CFD solution and iterations of the particle tracker.  NOTE:
  // the CFD time step must be an integral multiple of the particle
  // tracker time step.

  // ti: The current particle time step id.
  // ptc_int: The number of particle time steps between each CFD time step.
  // ptc_ftn: The number of particle time steps between each CFD results file.
  // user->ts: The number of CFD time steps between each results files.
  // time: the simulation time = it * ptc_dt
  // itime: The current particle time step id = it
  // time_o: time of the previously read CFD results file. = itime_o * ptc_dt
  // itime_o: particle time step id of the previously read CFD results file.
  // time_n: time of the current CFD results file. = itimen * ptc_dt
  // itime_n: particle time step id of the current CFD results file.
  // ptc_dt: The delta time of each particle time step.
  // user->dt: the delta time of each CFD time step.
  // user->tsdt: the delta time for ts timesteps.
  // IMPORTANT:
  // The avantages of using itime, itime_o and itime_n is that they are
  // integers so the same time ALWAY compare equal. This is not the case with
  // time, time_o and time_n which always need to tol equal compare.

  // Do setup things required for the restart.
  if (user->ptc_rstart > 0)
  {
    tistart = user->ptc_rstart * ptc_int;
    itime_n = tistart;
    time_n = tistart * ptc_dt;
    double searchRad = 3.0 * user->spatialIndexGrid->getMaxDelta();

    for (pn=0; pn<user->ptc_num; pn++)
    {
      if (Ptc[pn].stime <= GeomX::Tol(time_n))
      {
        Vector3d c(Ptc[pn].coor.x, Ptc[pn].coor.y, Ptc[pn].coor.z);
        auto result = user->rtreeGrid->getClosestEntity(c, &Ptc[pn].dist,
                                                                    searchRad);
        if (result.isValid())
        {
          Ptc[pn].geomLookupStatus = GEOM_SUCCESSFUL;
          Id(result).ijk(Ptc[pn].loc.i, Ptc[pn].loc.j, Ptc[pn].loc.k);
        }
        else
          Ptc[pn].geomLookupStatus = GEOM_FAILED;

        if (Ptc[pn].stime < GeomX::Tol(time_n))
          Ptc[pn].firstVel = false;
      }
    }
  }
  else
  {
    tistart = user->tis * ptc_int;

    if (user->setPartVel == 1)
    {
      int bi = 0;
      int fileidx = (tistart + ptc_ftn) / ptc_int;
      itime_o = tistart;
      time_o = tistart * ptc_dt;
      itime_n = fileidx * ptc_int;
      time_n = fileidx * ptc_dt * ptc_int;
      PetscPrintf(PETSC_COMM_WORLD,
                  "Reading file %d for time %f\n", fileidx, time_n);
      int ierr = VelocityReadin(user, fileidx, bi);
      if (ierr)
        return ierr;

      for (pn=0; pn<user->ptc_num; pn++)
      {
        if (Ptc[pn].stime <= GeomX::Tol(time_n))
        {
          if (Ptc[pn].geomLookupStatus == GEOM_NOT_DONE)
          {
            Ptc[pn].Vp.set(0, 1.0);
            ParticleLocation(user.get(), &Ptc[pn], 5, pn);
            Ptc[pn].Vp.set(0, 0.0);
          }
        }

        // The option has been set to calculate the initial velocity
        // for the particle to the fluid velocity plus the terminal velocity.
        if (Ptc[pn].stime == GeomX::Tol(time_o))
        {
          Location sloc;
          ISPointInCell(Ptc[pn].coor, user.get(), sloc, &Ptc[pn].coef, &Ptc[pn], pn);

          // Find the flow field velocity at the particle's location.
          VelocityInterpolation(user.get(), &Ptc[pn], time_o, time_o, time_n);

          // Add the terminal velocity.
          auto Vt = calcTermVel(user->DAccGrav, user->velScale, user, &Ptc[pn]);
          Vector3d const Vf = Ptc[pn].Vf1.get();
          Vector3d Vp(Vf + addTermVel(Vt, user->DAccGrav, Vector3d(0.0)));
          Ptc[pn].Vp.set(Vp);
        }
      }

      // Need to reread the first ufield file since each call to VelocityReadin
      // the current ufield to old and than reads in the new.
      ierr = VelocityReadin(user, user->tis, bi);
      if (ierr)
        return ierr;

      itime_n = tistart;
      time_n = tistart * ptc_dt;
    }

    ParticleTrajectoryIO2(user, Ptc, user->tis, user->tis * user->tsdt);
    ParticleTrajectoryIO3(user, Ptc, user->tis, user->tis * user->tsdt);
  }



  // DEBUG START
#if 0
  Location const DBloc(Ptc[0].loc);
  Cmpnts const DBcoor(Ptc[0].coor);

  {
    auto fd = fopen("FV.profile.dat", "w");
    fprintf(fd, "ti\ttime\tFvx\tFvy\tFvz\n");
    fclose(fd);
  }
#endif
  // DEBUG END


  for (ti=tistart+1; ti<=user->tie*ptc_int; ti++) 
  {
    // If the file named "STOP" appears in the directory, immediately
    // exit the program.
    auto fd = fopen("STOP", "r");
    if (fd)
    {
      PetscPrintf(PETSC_COMM_WORLD,
                  "STOP file exist, program is stopping\n");
      fclose(fd);
      return false;
    }

    itime = ti;
    time = ti * ptc_dt;

    if (ti % ptc_ftn == 1)
    {
      int fileidx = (ti - 1 + ptc_ftn) / ptc_int;
      itime_o = itime_n;
      time_o = time_n;
      itime_n = fileidx * ptc_int;
      time_n = fileidx * ptc_dt * ptc_int;
      
      PetscPrintf(PETSC_COMM_WORLD,
                  "Reading file %d for time %f\n", fileidx, time);
        
      // for (bi=0; bi < blocknumber; bi++)
      {
        bi = 0;
        int ierr =
          // VelocityReadin(&user[bi], fileidx, bi);
          VelocityReadin(user, fileidx, bi);
        if (ierr)
          return ierr;
      }

      if (user->useDebugDrag)
        user->debugDrag->swapFile(fileidx);
    }

    PetscPrintf(PETSC_COMM_WORLD, "Particle Tracking time %i - %f\n", ti, time);



    // DEBUG START
#if 0
    // uses location, coor, sets coef
    Ptc[0].loc = DBloc;
    Ptc[0].coor = DBcoor;
    Location sloc;
    ISPointInCell(Ptc[0].coor, user, sloc, &Ptc[0].coef, &Ptc[0], 0);

    // uses coef, loc, sets P->Vf
    // find the particle velocity from flow field
    VelocityInterpolation(user, &Ptc[0], time, time_o, time_n);

    {
      auto fd = fopen("FV.profile.dat", "a");
      fprintf(fd, "%d\t%6.3f\t%15e\t%15e\t%15e\n", ti, time, Ptc[0].Vf.x, Ptc[0].Vf.z, Ptc[0].Vf.z);
      fclose(fd);
    }
#endif
    // DEBUG END



    // Ali's covid breathing simulation - Check if particles need to be added.
    // Not using this routine now since it generates random locations.  In
    // order to get repeatable results, I generated a file of particle
    // definitions that used the same random process so they will be the
    // same for any simulation.
    // if (user->useBreath == 1)
    //   breathAddParticles(ti, time, ptc_dt, user, Ptc);
    
    for (pn=0; pn<user->ptc_num; pn++) 
    {
      if (Ptc[pn].stime > GeomX::Tol(time) || Ptc[pn].partStat != INSIDE_ACTIVE)
        continue;

      // For particles that are starting at this time, they don't need to
      // have their trajectory calculated since their position at their
      // start time is their current position.

      // However, all particles need to have an i,j,k value assigned to them
      // before any real processing can continue so find them for all whose
      // time has not yet come and for whose who are starting with this
      // interation.
      if (Ptc[pn].stime <= GeomX::Tol(time))
      {
        if (Ptc[pn].geomLookupStatus == GEOM_NOT_DONE)
          ParticleLocation(user.get(), &Ptc[pn], 5, pn);

        // The option has been set to calculate the initial velocity
        // for the particle to the fluid velocity plus the terminal velocity.
        if (user->setPartVel == 1 && Ptc[pn].stime == GeomX::Tol(time))
        {
          Location sloc;
          ISPointInCell(Ptc[pn].coor, user.get(), sloc, &Ptc[pn].coef, &Ptc[pn], pn);

          // Find the flow field velocity at the particle's location.
          VelocityInterpolation(user.get(), &Ptc[pn], time, time_o, time_n);

          // Add the terminal velocity.
          auto Vt = calcTermVel(user->DAccGrav, user->velScale, user, &Ptc[pn]);
          Vector3d Vf = Ptc[pn].Vf1.get();
          Vector3d Vp = Vf + addTermVel(Vt, user->DAccGrav, Vector3d(0.0));
          Ptc[pn].Vp.set(Vp);
        }
      }
    
      // Ali's covid breathing simulation
      // Calculate new diameters for particles
      if (user->useBreath == 1)
      {
        breathCalcPartDia(time, Ptc[pn]);
        if (Ptc[pn].DDiaPt <= GeomX::Tol(0.0))
        {
          Ptc[pn].DDiaPt = 0;
          updatePartStatus(pn, ZERO_DIAMETER, Ptc[pn].loc,
                           Vector3d(Ptc[pn].coef.x,
                                    Ptc[pn].coef.y, Ptc[pn].coef.z), &Ptc[pn]);
          continue;
        }
      }


      Location sloc;
      ISPointInCell(Ptc[pn].coor, user.get(), sloc, &Ptc[pn].coef, &Ptc[pn], pn);

      // Find the flow field velocity at the particle's location.
      VelocityInterpolation(user.get(), &Ptc[pn], time, time_o, time_n);
    
      ParticleTrajectory(user, &Ptc[pn], ti, pn, time, time_n, time_o);

      ParticleLocation(user.get(), &Ptc[pn], 5, pn);
    }


    if (user->output_ts == 0)
    {
      if (fmod((double)ti / ptc_int, user->ts) == 0.0)
        ParticleTrajectoryIO2(user, Ptc, ti / ptc_int, time);
    }
    else if (fmod((double)ti / ptc_int, user->output_ts) == 0.0)
      ParticleTrajectoryIO2(user, Ptc, ti / ptc_int, time);

    if (user->useDebugTraj > 0 || itime == itime_n)
      ParticleTrajectoryIO3(user, Ptc, itime_n / ptc_int, time);
  }



  
  DAVecRestoreArray(user->fda, Coor, &user->coor);

  DAVecRestoreArray(user->da,  user->Rad,  &user->rad);
  DAVecRestoreArray(user->fda, user->Cent, &user->cent);

  DAVecRestoreArray(user->fda, user->Ucat,   &user->ucat);
  DAVecRestoreArray(user->fda, user->Ucat_o, &user->ucat_o);

  return 0;
}

PetscTruth ISPointInCell(Cmpnts p, UserCtx * user,
                         Location loc, Cmpnts *coef, Particle *P, int pId,
                         bool printErr)
{
  // Don't process particles unless they are active.
  if (P->partStat != INSIDE_ACTIVE)
    return PETSC_FALSE;
  
  PetscInt const i = P->loc.i;
  PetscInt const j = P->loc.j;
  PetscInt const k = P->loc.k;


  // For my code (WRO), this routine calculates P->coef ONLY.

  // You stupid ass! You know you shouldn't change working code unless
  // you absolutely need to.  Remember this lesson and stop it!  WRO 2021-03-09
  
#if 0
  double dis[6];
  Vector3d pv(p.x, p.y, p.z);
  std::vector<Vector3d> fcoor(4);

  
#if 0
  int vId = 0;
  Vector3d vcoor;
  PetscPrintf(PETSC_COMM_WORLD, "\n");
  for (int ii=0; ii<8; ++ii)
  {
    user->spatialIndexGrid->getVertCoords(ii, i, j, k, user->coor, vcoor);
    PetscPrintf(PETSC_COMM_WORLD, "v,%d,%f,%f,%f\n",
                ++vId, vcoor(0), vcoor(1), vcoor(2));
  }
  PetscPrintf(PETSC_COMM_WORLD, "v,%d,%f,%f,%f\n", ++vId, pv(0), pv(1), pv(2));
#endif  

  // Find coef in x which we assume to be in the k direction.
  // Uses faces 5 and 4
  user->spatialIndexGrid->getFaceCoords(5, i, j, k, user->coor, fcoor);
  dis[5] = sqrt(Geom::facePointDistanceSquared(fcoor, pv));

  user->spatialIndexGrid->getFaceCoords(4, i, j, k, user->coor, fcoor);
  dis[4] = sqrt(Geom::facePointDistanceSquared(fcoor, pv));

  // Find coef in y which we assume to be in the i direction.
  // Uses faces 1 and 0
  user->spatialIndexGrid->getFaceCoords(1, i, j, k, user->coor, fcoor);
  dis[1] = sqrt(Geom::facePointDistanceSquared(fcoor, pv));

  user->spatialIndexGrid->getFaceCoords(0, i, j, k, user->coor, fcoor);
  dis[0] = sqrt(Geom::facePointDistanceSquared(fcoor, pv));

  // Find coef in z which we assume to be in the j direction.
  // Uses faces 3 and 2
  user->spatialIndexGrid->getFaceCoords(3, i, j, k, user->coor, fcoor);
  dis[3] = sqrt(Geom::facePointDistanceSquared(fcoor, pv));

  user->spatialIndexGrid->getFaceCoords(2, i, j, k, user->coor, fcoor);
  dis[2] = sqrt(Geom::facePointDistanceSquared(fcoor, pv));

  coef->x = dis[5] / (dis[5] + dis[4]);
  coef->y = dis[1] / (dis[1] + dis[0]);
  coef->z = dis[3] / (dis[3] + dis[2]);

//PetscPrintf(PETSC_COMM_WORLD, "coef #1: %f,%f,%f\n", coef->x,coef->y,coef->z);

  return PETSC_TRUE;
  
#else

  Cmpnts cell_v[8];

  cell_v[0] = user->coor[k][j][i];
  cell_v[1] = user->coor[k][j][i+1];
  cell_v[2] = user->coor[k][j+1][i+1];
  cell_v[3] = user->coor[k][j+1][i];

  cell_v[4] = user->coor[k+1][j][i];
  cell_v[5] = user->coor[k+1][j][i+1];
  cell_v[6] = user->coor[k+1][j+1][i+1];
  cell_v[7] = user->coor[k+1][j+1][i];


  PetscReal dis[6];
  // My face 0
  PFaceDis(p, cell_v[1], cell_v[2], cell_v[6], cell_v[5], &dis[0]);
  int lastCmp = 0;
  if (dis[0] <= GeomX::Tol10<double,10>(0.0)) {
    // My face 1
    PFaceDis(p, cell_v[0], cell_v[3], cell_v[7], cell_v[4], &dis[1]);
    lastCmp = 1;
    if (dis[1] >= GeomX::Tol10<double,10>(0.0)) {
      // coef->x = dis[1] / (dis[1] - dis[0]);  original direction
      coef->y = dis[1] / (dis[1] - dis[0]);
      if (coef->y < 0.0)
        coef->y = 0.0;
      if (coef->y > 1.0)
        coef->y = 1.0;
      // My face 2
      PFaceDis(p, cell_v[3], cell_v[7], cell_v[6], cell_v[2], &dis[2]);
      lastCmp = 2;
      if (dis[2] <= GeomX::Tol10<double,10>(0.0)) {
        // My face 3
        PFaceDis(p, cell_v[0], cell_v[4], cell_v[5], cell_v[1], &dis[3]);
        lastCmp = 3;
        if (dis[3] >= GeomX::Tol10<double,10>(0.0)) {
	  // coef->y = dis[3] / (dis[3] - dis[2]);  original direction
	  coef->z = dis[3] / (dis[3] - dis[2]);
          if (coef->z < 0.0)
            coef->z = 0.0;
          if (coef->z > 1.0)
            coef->z = 1.0;
          // My face 4
          PFaceDis(p, cell_v[4], cell_v[5], cell_v[6], cell_v[7], &dis[4]);
          lastCmp = 4;
          if (dis[4] <= GeomX::Tol10<double,10>(0.0)) {
            // My face 5
            PFaceDis(p, cell_v[0], cell_v[1], cell_v[2], cell_v[3], &dis[5]);
            lastCmp = 5;
            if (dis[5] >= GeomX::Tol10<double,10>(0.0)) {
	      // coef->z = dis[5] / (dis[5] - dis[4]);  original direction
	      coef->x = dis[5] / (dis[5] - dis[4]);
              if (coef->x < 0.0)
                coef->x = 0.0;
              if (coef->x > 1.0)
                coef->x = 1.0;
              // PetscPrintf(PETSC_COMM_WORLD, "coef #2: %f,%f,%f\n", coef->x,coef->y,coef->z);
	      return PETSC_TRUE;
	    }
	  }
	}
      }
    }
  }

  if (printErr)
  {
    user->globalError |= COEF_ERROR;
    PetscPrintf(PETSC_COMM_WORLD,
                "Coef calc error:  pId: %d  i,j,k: %d,%d,%d\n", pId, i, j, k);

    Vector3d pv(p.x, p.y, p.z);
    PetscPrintf(PETSC_COMM_WORLD, "v,%d,%.15g,%.15g,%.15g\n", 10, pv(0), pv(1), pv(2));

    static char const * cmpType[] = { "<=", "=>", "<=", ">=", "<=", ">=" };

    PetscPrintf(PETSC_COMM_WORLD, "lastCmp: %d\n", lastCmp);
    for (int i=0; i<=lastCmp; ++i)
      PetscPrintf(PETSC_COMM_WORLD, "dis[%d](%s): %.15g\n", i, cmpType[i], dis[i]);


    int vId = 0;
    Vector3d vcoor;
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    for (int ii=0; ii<8; ++ii)
    {
      user->spatialIndexGrid->getVertCoords(ii, i, j, k, user->coor, vcoor);
      PetscPrintf(PETSC_COMM_WORLD, "v,%d,%.15g,%.15g,%.15g\n",
                  ++vId, vcoor(0), vcoor(1), vcoor(2));
    }
  }
  
  return PETSC_FALSE;

#endif
  
}


template < typename T >
void VecAB(T & c, T const & a, T const & b)
{
  c.x = b.x - a.x;
  c.y = b.y - a.y;
  c.z = b.z - a.z;
}

template < typename T >
void CROSS(T & c, T const & a, T const & b)
{
  c.x = a.y * b.z - a.z * b.y;
  c.y = a.z * b.x - a.x * b.z;
  c.z = a.x * b.y - a.y * b.x;
}

template < typename T >
double DOT(T const & a, T const & b)
{
  return a.x*b.x + a.y*b.y + a.z*b.z;
}


PetscErrorCode PFaceDis(Cmpnts p, Cmpnts p0, Cmpnts p1, Cmpnts p2, Cmpnts p3, PetscReal *dis)
{
  Cmpnts p02, p13;

  Cmpnts pc;
  Cmpnts norm;
  PetscReal length;

  VecAB(p02, p0, p2);
  VecAB(p13, p1, p3);

  CROSS(norm, p02, p13);

  length = sqrt(norm.x * norm.x + norm.y * norm.y + norm.z * norm.z);

  pc.x = 0.25 * (p0.x + p1.x + p2.x + p3.x);
  pc.y = 0.25 * (p0.y + p1.y + p2.y + p3.y);
  pc.z = 0.25 * (p0.z + p1.z + p2.z + p3.z);

  norm.x /= length; norm.y /= length; norm.z /= length;

  Cmpnts pc_p;
  VecAB(pc_p, pc, p);
  *dis = DOT(pc_p, norm);
  *dis = (*dis == GeomX::Tol(0.0)) ? 0.0 : *dis;
  return 0;
}

PetscErrorCode ParticleLocation(UserCtx * user, Particle *P,
                                PetscInt iradial, int pId, bool outsideIsAnError)
{
 if (P->partStat != INSIDE_ACTIVE)
    return true;

#if 1
  Vector3d c;
  int i=0, j=0, k=0;
  Geom::Id result;
  
  if (P->geomLookupStatus == GEOM_SUCCESSFUL)
  {
      i = P->loc.i;
      j = P->loc.j;
      k = P->loc.k;
  }

  if (P->geomLookupStatus == GEOM_NOT_DONE)
  {
    if (wroDebugPID == -1 || wroDebugPID == pId)
      PrintDebugInfo01(true, P, pId, 0, Vector3d(0.0), Vector3d(0.0),
                       Vector3d(0.0),
                       Vector3d(0.0),
                       Vector3d(P->coor.x, P->coor.y, P->coor.z));

    c = Vector3d(P->coor.x, P->coor.y, P->coor.z);
    double searchRad = 3.0 * user->spatialIndexGrid->getMaxDelta();
    result = user->rtreeGrid->getClosestEntity(c, &P->dist, searchRad);

    if (result.isValid())
    {
      P->geomLookupStatus = GEOM_SUCCESSFUL;
      Id(result).ijk(i, j, k);
      P->loc.i = i;
      P->loc.j = j;
      P->loc.k = k;
    }
    else
      P->geomLookupStatus = GEOM_FAILED;
  }

  if (P->geomLookupStatus == GEOM_FAILED)
  {
    c = Vector3d(P->coor.x, P->coor.y, P->coor.z);
    result = user->rtreeGrid->getClosestEntity(c, &P->dist);

    if (result.isValid())
    {
      P->geomLookupStatus = GEOM_SUCCESSFUL;
      Id(result).ijk(i, j, k);
      P->loc.i = i;
      P->loc.j = j;
      P->loc.k = k;
    }
    else
      P->geomLookupStatus = GEOM_FAILED;
  }

  if (P->geomLookupStatus == GEOM_FAILED || (outsideIsAnError && P->dist != 0.0))
  {
    updatePartStatus(pId, LOOKUP_FAILED, P->loc,
                     Vector3d(P->coef.x, P->coef.y, P->coef.z), P);
    P->coef.x = P->coef.y = P->coef.z = 0.;
    return true;
  }

  if (P->dist != 0.0)
  {
    P->partStat = LOOKUP_FAILED;
    return OUTSIDE_UNKNOWN;
  }


  // DEBUGGING - is the timestep too big.
  {
    int maxDelta = std::max(abs(i - P->loc.i), abs(j - P->loc.j));
        maxDelta = std::max(maxDelta,          abs(k - P->loc.k));
    if (maxDelta > 1)
      PetscPrintf(PETSC_COMM_WORLD,
                  "pId: %d  Timestep too big!  delta: %d\n", pId, maxDelta);
  }

    
  // if (!P->firstVel && vmag < 0.001)
  if (!P->firstVel && P->Vp.getv().mag() < 1.0e-6)
  {
    Location loc;
    loc.i = i;
    loc.j = j;
    loc.k = k;
    auto const & c = P->coor;
    updatePartStatus(pId, ZERO_VEL_LOC, loc, Vector3d(c.x, c.y, c.z), P);
    P->coef.x = P->coef.y = P->coef.z = 0.;
    return true;
  }
  
  P->firstVel = false;
  
  return false;
  
#else

  Cmpnts p = P->coor;
  Location loc = P->loc;

  DALocalInfo	info = user->info;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt lis, lie, ljs, lje, lks, lke;

  PetscInt i=0, j=0, k=0;

  Cmpnts ***cent;
  DAVecGetArray(user->fda, user->Cent, &cent);

  lis = PetscMax(0, loc.i - iradial);
  lie = PetscMin(mx-2, loc.i + iradial);
  ljs = PetscMax(0, loc.j - iradial);
  lje = PetscMin(my-2, loc.j + iradial);
  lks = PetscMax(0, loc.k - iradial);
  lke = PetscMin(mz-2, loc.k + iradial*2);

  PetscTruth ISInside = PETSC_FALSE;
  Location sloc;
  Cmpnts v;
  for (k=lks; k<lke && !ISInside; k++) {
    for (j=ljs; j<lje && !ISInside; j++) {
      for (i=lis; i<lie && !ISInside; i++) {
	VecAB(v, p, cent[k][j][i]);
	if (Length(v) < user->rad[k][j][i]) {
	  sloc.i = i; sloc.j = j; sloc.k = k;

	  ISInside = ISPointInCell(p, user, sloc, &P->coef, P, pId);
	  if (ISInside) {
	    P->loc.i = i; P->loc.j = j; P->loc.k = k;
	  }
	}
      }
    }
  }

  if (j > user->JM-3 ) ISInside = PETSC_FALSE; /* mli */
  if (!ISInside) {
    PetscPrintf(PETSC_COMM_WORLD, "Outside domain!\n");
    P->coef.x = 0; P->coef.y = 0; P->coef.z = 0.;
  }
  
  DAVecRestoreArray(user->fda, user->Cent, &cent);
  return ISInside;

#endif

}

#define TriLinear(v000, v001, v010, v100, v011, v110, v101, v111, coef) \
        v000 * (1.0-coef.x) * (1.0-coef.z) * (1.0-coef.y) + \
        v100 *      coef.x  * (1.0-coef.z) * (1.0-coef.y) + \
        v010 * (1.0-coef.x) *      coef.z  * (1.0-coef.y) + \
        v001 * (1.0-coef.x) * (1.0-coef.z) *      coef.y  + \
        v101 *      coef.x  * (1.0-coef.z) *      coef.y  + \
        v011 * (1.0-coef.x) *      coef.z  *      coef.y  + \
        v110 *      coef.x  *      coef.z  * (1.0-coef.y) + \
        v111 *      coef.x  *      coef.z  *      coef.y

PetscErrorCode VelocityInterpolation(UserCtx * user,
                                     Particle *P, PetscReal time,
                                     PetscReal time_o, PetscReal time_n, bool VF0update)
{
  if (P->partStat != INSIDE_ACTIVE)
    return 0;

  // Move current to previous since this routine calculates current.
  if (VF0update)
    P->Vf0.set(P->Vf1.get());
  
  PetscReal interpolation;
  Cmpnts coef = P->coef;
  interpolation = (time - time_o) / (time_n - time_o);
  
  Cmpnts V000, V001, V010, V100, V011, V110, V101, V111;

  auto const * const * const * const ucat = user->ucat;
  auto const * const * const * const ucat_o = user->ucat_o;
  Vector3ds vold, vnew;

  // The velocities at the extreme of the grid are set to zero.  I'm going
  // to setup indices to bypass those problems.
  PetscInt i, j, k, i1, j1, k1;

  i = P->loc.i;
  j = P->loc.j;
  k = P->loc.k;
  i1 = i + 1;
  j1 = j + 1;
  k1 = k + 1;

  if (i == 0)
    i = i1;
  if (i1 == user->IM-1)
    i1 = i;

  if (j == 0)
    j = j1;
  if (j1 == user->JM-1)
    j1 = j;

  if (k == 0)
    k = k1;
  if (k1 == user->KM-1)
    k1 = k;

/*   DAVecGetArray(user->fda, user->Ucat, &ucat); */

  V000 = ucat[k][j][i];
  V001 = ucat[k][j][i1];
  V010 = ucat[k][j1][i];
  V100 = ucat[k1][j][i];
  V011 = ucat[k][j1][i1];
  V110 = ucat[k1][j1][i];
  V101 = ucat[k1][j][i1];
  V111 = ucat[k1][j1][i1];

/*   DAVecRestoreArray(user->fda, user->Ucat, &ucat); */
  vnew.set(TriLinear(V000.x, V001.x, V010.x, V100.x,
		     V011.x, V110.x, V101.x, V111.x, coef),

           TriLinear(V000.y, V001.y, V010.y, V100.y,
		     V011.y, V110.y, V101.y, V111.y, coef),

           TriLinear(V000.z, V001.z, V010.z, V100.z,
		     V011.z, V110.z, V101.z, V111.z, coef));


/*   DAVecGetArray(user->fda, user->Ucat_o, &ucat); */

  V000 = ucat_o[k][j][i];
  V001 = ucat_o[k][j][i1];
  V010 = ucat_o[k][j1][i];
  V100 = ucat_o[k1][j][i];
  V011 = ucat_o[k][j1][i1];
  V110 = ucat_o[k1][j1][i];
  V101 = ucat_o[k1][j][i1];
  V111 = ucat_o[k1][j1][i1];

/*   DAVecRestoreArray(user->fda, user->Ucat_o, &ucat); */
  vold.set(TriLinear(V000.x, V001.x, V010.x, V100.x,
		     V011.x, V110.x, V101.x, V111.x, coef),

           TriLinear(V000.y, V001.y, V010.y, V100.y,
		     V011.y, V110.y, V101.y, V111.y, coef),

           TriLinear(V000.z, V001.z, V010.z, V100.z,
		     V011.z, V110.z, V101.z, V111.z, coef));

  P->Vf1.set(vold.getv() + interpolation * (vnew.getv() - vold.getv()));

  return 0;
}

PetscErrorCode ParticleTrajectory(std::shared_ptr< UserCtx > & user,
                                  Particle * const P, PetscInt ti, int pId,
                                  double time, double time_o, double time_n)
{
  if (P->partStat != INSIDE_ACTIVE)
    return 0;

  double searchRad = 3.0 * user->spatialIndexGrid->getMaxDelta();

  // This variable is used thoughout this subroutine.  Make sure to only
  // use is locally.

  Geom::Id result;
  
  P->etime = time;

  // Nondimensional scaling factors.  If you are running a dimensional
  // simulation, these variables must to 1 which is the default.
  double const L = user->lengScale;
  double const U = user->velScale;
  double const T = L / U;
  
  // If the current run is nondimensional, the variables must be changed
  // to dimensioned values to get the correct solution in this routine.
  // All dimensianed variables start with D.

  Vector3ds const P0(P->coor);
  Location const loc0(P->loc);
  // Fluid velocity
  Vector3ds const Vf0(P->Vf0);
  Vector3ds const Vf1(P->Vf1);
  // Particle velocity
  Vector3ds const Vp0(P->Vp);
  // Vpf0 is the particle's velocity w.r.t. the fluid.
  Vector3ds const Vpf0(Vp0.getv() - Vf0.getv());

  Vector3d Vpf1(0.0);
  Vector3d Vp1(0.0);
  Vector3d P1(0.0);
  Location loc1;
  double DMasPt = 0.0;
  Vector3d DAccGrav(user->DAccGrav);
  double rey = 0.0;
  P->Rey = 0.0;
  P->RespTime = 0.0;
  double Cd = 0.0;
  Vector3d DFdrag1(0.0);
  Vector3d DFdrag2(0.0);
  Vector3d DFdrag3(0.0);
  Vector3d DFdrag(0.0);
  Vector3d DFtot(0.0);
  Vector3d DAcc(0.0);
  double Vt = 0.0;




  // Get results for RK method to compare with current calculations.
  RKY RKresult;
  if (user->calcMethod == 4 || user->RKtest)
  {
    RKparms rkparms;
    rkparms.u = user.get();
    rkparms.p = P;
    rkparms.pId = pId;
    rkparms.time = time;
    rkparms.time_o = time_o;
    rkparms.time_n = time_n;

    std::unique_ptr<Particle> saveP;
    if (user->RKtest)
    {
      saveP.reset(new Particle);
      memcpy(saveP.get(), P, sizeof(Particle));
    }

    RKresult = RungeKutta(rkparms, ptc_dt, time, RKY(P->Vp, P->coor), calcParticleAcceleration);
    if (RKresult.getOutside())
      int debug = 1;

#if 0
    // I was playing with adaptive step size here. Didn't see much improvement so I didn't
    // implement any.
    RKresult1 = RungeKuttaAdaptive(rkparms, ptc_dt, time, RKY(P->Vp, P->coor), calcParticleAcceleration);
    memcpy(P, saveP.get(), sizeof(Particle));
    P->geomLookupStatus = GEOM_NOT_DONE;

    RKresult2 = RungeKuttaAdaptive(rkparms, ptc_dt/2.0, time, RKY(P->Vp, P->coor),
                                   calcParticleAcceleration);
    RKresult3 = RungeKuttaAdaptive(rkparms, ptc_dt/2.0, time+ptc_dt/2.0, RKresult2,
                                   calcParticleAcceleration);
    RKerr = RKresult3 - RKresult;
    memcpy(P, saveP.get(), sizeof(Particle));
    P->geomLookupStatus = GEOM_NOT_DONE;

    RKresult4 = RungeKuttaAdaptive(rkparms, ptc_dt/4.0, time,                RKY(P->Vp, P->coor),
                                   calcParticleAcceleration);
    RKresult5 = RungeKuttaAdaptive(rkparms, ptc_dt/4.0, time+ptc_dt/4.0,     RKresult4,
                                   calcParticleAcceleration);
    RKresult6 = RungeKuttaAdaptive(rkparms, ptc_dt/4.0, time/2.0,            RKresult5,
                                   calcParticleAcceleration);
    RKresult7 = RungeKuttaAdaptive(rkparms, ptc_dt/4.0, time+3.0*ptc_dt/4.0, RKresult6,
                                   calcParticleAcceleration);
    RKerr = RKresult7 - RKresult;
    RKerr2 = RKresult7 - RKresult3;
    memcpy(P, saveP.get(), sizeof(Particle));
    P->geomLookupStatus = GEOM_NOT_DONE;
#endif

    if (user->RKtest)
    {
      P->geomLookupStatus = GEOM_NOT_DONE;
      memcpy(P, saveP.get(), sizeof(Particle));
    }

    else
    {
      Vp1 = RKresult.getv();
      P1 = RKresult.getx();
    }
  }

  if (user->calcMethod != 0 || (user->calcMethod != 0 && user->RKtest))
  {
    int saveCalcMethod;
    if (user->RKtest)
    {
      saveCalcMethod = user->calcMethod;
      user->calcMethod = 1;
    }
      
    // Calculate the particle mass.
    double DVolPt = Geom::PI * GeomX::pow(P->DDiaPt, 3) / 6;
    DMasPt = P->DDenPt * DVolPt;

    // Calculate the drag force
    P->Rey = rey = user->DDenFd * P->DDiaPt * (U*Vpf0.getv().mag()) / user->DVisFd;
    P->RespTime = (P->DDenPt * P->DDiaPt * P->DDiaPt) / (18.0 * user->DVisFd);

    if (user->calcMethod == 3)
    {
      DFdrag3 = dragKineticGas(user->DVisFd, P->DDiaPt, U*Vpf0.getv());

      // Written this way to work with debugging code.
      // Used for debugging only.
      DFdrag = DFdrag3;
      DFtot = DFdrag + DMasPt * DAccGrav;

      // Calculate the acceleration and velocity that would be due to that a.
      // then calculate the new particle velocity w.r.t. the fluid.
      DAcc = DFdrag / DMasPt + DAccGrav;
      Vpf1 = DAcc * (T*ptc_dt) / U + Vpf0;
      Vp1 = Vf1.getv() + Vpf1;

      P1 = ((0.5 * DAcc * (T*ptc_dt)) + (U*Vp1)) * (T*ptc_dt) / L + P0;
    }
    else
    {
      // At very small length scales, this method becomes unstable because we
      // are in the machine zero regime so we are basically calculating only
      // with rounding error.  In this case, set Vpf1 to 0 which may be done
      // indirectly by setting the drag to zero. The number I choose as the
      // break point is based on a breath/noface looking at why one particle
      // didn't work correctly. The time stamp was 27.14 s.
      if (rey <= 5.0e-6)
      {
        DFdrag1 = Vector3d(0.0);
        user->globalError |= FORCE_REL_VEC_TO_ZERO;
        // PetscPrintf(PETSC_COMM_WORLD,
        //         "Particles relative velocity set to zero: pId: %d\n", pId);
      }

      else if (rey <= 0.001) // This stops divide by zero and gets a
                             // valid result.
        DFdrag1 = -3.0 * Geom::PI * user->DVisFd * P->DDiaPt * (U*Vpf0.getv());

      else
      {
        Cd = ReysToCd::conv(rey);
        DFdrag1 = -0.125 * user->DDenFd * Cd * Geom::PI
          * GeomX::pow(P->DDiaPt, 2) * (U*Vpf0.getv().mag()) * (U*Vpf0.getv());
      }

      DFdrag2 = dragKineticGas(user->DVisFd, P->DDiaPt, U*Vpf0.getv());

      if (user->calcMethod == 1)
      {
        // PetscPrintf(PETSC_COMM_WORLD, "pId: %d  DFdrag1.mag: %e  DFdrag2.mag: %e  ratio: %e\n",
        //             pId, DFdrag1.mag(), DFdrag2.mag(), DFdrag2.mag() / DFdrag1.mag());
        if (DFdrag1.mag() < DFdrag2.mag())
          DFdrag = DFdrag1;
        else
          DFdrag = DFdrag2;
      }
      else
        DFdrag = DFdrag2;
    
      // Used for debugging only.
      DFtot = DFdrag + DMasPt * DAccGrav;

      // Calculate the acceleration and velocity that would be due to that a.
      // then calculate the new particle velocity w.r.t. the fluid.
      DAcc = DFdrag / DMasPt + DAccGrav;
      Vpf1 = DAcc * (T*ptc_dt) / U + Vpf0;
      Vp1 = Vf1.getv() + Vpf1;

      // Calculate the terminal velocity.
      Vt = calcTermVel(DAcc, U, user, P);

      // If the new velocity is > the terminal velocity, modify Vpf1 to have
      // the magnitude of Vt.  P1 is the new particle coordinate.
      if (Vpf1.mag() > Vt)
      {
        // if (wroDebugPID == -2)
        //   PetscPrintf(PETSC_COMM_WORLD, "Part Id: %d  Speed limit exceded: "
        //                                 "Dia: %f  Re: %f  Cd: %f  Vpf0: %f  "
        //                                 "Vpf1: %f  Vt: %f\n",
        //    pId, P->dia*1e6, rey, Cd, (U*Vpf0.mag()), (U*Vpf1.mag()), (U*Vt));

        Vpf1 *= Vt / Vpf1.mag();
        Vp1 = Vf1.getv() + Vpf1;
        P1 = (U*Vp1) * (T*ptc_dt) / L + P0;
      }
      else
      {
        P1 = ((0.5 * DAcc * (T*ptc_dt)) + (U*Vp1)) * (T*ptc_dt) / L + P0;
      }
    }

    
    if (user->useDebugDrag)
    {
      double dist;
      Location loc1;
      result = user->rtreeGrid->getClosestEntity(P1, &dist, searchRad);
      if (result.isValid())
        Id(result).ijk(loc1.i, loc1.j, loc1.k);

      user->debugDrag->add(user, pId, time, Vf0.get(), Vp0.get(), Vpf0, P->DDiaPt,
                           user->DDenFd, P->DDenPt, DMasPt, DAccGrav(2),
                           rey, P->RespTime, Cd, DFdrag1, DFdrag2, DFdrag3,
                           DFdrag, DFtot, DAcc, DAcc * (T*ptc_dt) / U + Vpf0,
                           Vt, Vf1.get(), Vp1, Vpf1, P0, loc0, P1, loc1, P->partStat);
    }
    
    if (wroDebugPID == -1 || wroDebugPID == pId)
    {
      PrintDebugInfo01(false, P, pId, ti, DAcc, Vf1.get(), Vpf1, Vp1, P1);

      auto Pdel = P1 - P0;
      PrintDebugInfo02(pId, P0, ptc_dt, DVolPt, DMasPt, DAccGrav, Vf1, Vp0.get(),
                       Vpf0, rey, Cd, DFdrag, DFtot, DAcc, DAcc * (T*ptc_dt),
                       Vpf1, Vt, Vp1, Pdel, P1);
    }

    if (user->RKtest)
      user->calcMethod = saveCalcMethod;
  }

  // Eulerian/Lagrangian particle tracking method.
  if (user->calcMethod == 0 || (user->calcMethod == 0 && user->RKtest))
  {
    // Update velocity
    Vp1 = Vf1;

    // Calculate settling velocity
    if (DAccGrav.mag() != GeomX::Tol(0,0))
    {
      Vt = calcTermVel(DAccGrav, U, user, P);
      Vp1 = Vf1.getv() + addTermVel(Vt, DAccGrav, Vector3d(0.0));
    }

    // Update position
    P1 = P0.getv() + (U*Vp1) * (T*ptc_dt) / L;


    if (wroDebugPID == -1 || wroDebugPID == pId)
    {
      auto Vpf1 = Vpf0.getv() + Vp1;
      PrintDebugInfo01(false, P, pId, ti, Vector3d(0.0), Vf1, Vpf1, Vp1, P1);
    }
  }

  
  // See if particle is outside the grid.
  P->geomLookupStatus = GEOM_FAILED;
  raytrace2(P0, P1, pId, P, user, result, loc1, searchRad);

  
  if (user->useIBs && P->partStat == INSIDE_ACTIVE)
  {
    // If the new location is inside an ib, move the particle onto the
    // ib surface where it passed through.
    double dist;
    result = user->rtreeSurfs->getClosestEntity(P1, &dist, searchRad);
    if (result.isValid())
    {
      // If we are on the surface, mark particle as on an ib boundary.
      if (dist == GeomX::Tol(0.0))
      {
        auto result = user->rtreeGrid->getClosestEntity(P1, &dist, searchRad);
        if (result.isValid() && dist <= GeomX::Tol(0.0))
        {
          Id(result).ijk(loc1.i, loc1.j, loc1.k);
          updatePartStatus(pId, ON_IB_SURFACE, loc1, P1, P);
          P->geomLookupStatus = GEOM_SUCCESSFUL;
        }
      }

      // Don't allow the particle to move a further distance than it
      // is from P0.
      if ((P1-P0).mag() < dist)
      {
      }
      
      // See which side of the ib surface we are on.  This is possible since
      // the faces must be oriented with a right handed wrapping so that the
      // area normal is outward pointing.
      else
      {
        auto coors = user->surfs->getTriCoor(result);
        Vector3d area, centroid;
        Geom::triangleAreaCentroid(coors[0], coors[1], coors[2],
                                   area, centroid);
        area.normalize();
        auto VecDir = centroid - P1;
        VecDir.normalize();
        auto dir = area.dot(VecDir);

        if (wroDebugPID == -1 || wroDebugPID == pId)
        {
          PetscPrintf(PETSC_COMM_WORLD, "\nParticleTrajectory\n");
          PetscPrintf(PETSC_COMM_WORLD, "dist: %f\n", dist);
          PetscPrintf(PETSC_COMM_WORLD, "view,%f,%f,%f\n",
                      area(0), area(1), area(2));

          PetscPrintf(PETSC_COMM_WORLD, "v,1,%f,%f,%f\n",
                      coors[0](0), coors[0](1), coors[0](2));
          PetscPrintf(PETSC_COMM_WORLD, "v,2,%f,%f,%f\n",
                      coors[1](0), coors[1](1), coors[1](2));
          PetscPrintf(PETSC_COMM_WORLD, "v,3,%f,%f,%f\n",
                      coors[2](0), coors[2](1), coors[2](2));

          PetscPrintf(PETSC_COMM_WORLD, "v,4,%f,%f,%f\n",
                      centroid(0), centroid(1), centroid(2));
          auto tmp = centroid + area * searchRad / 5.0;
          PetscPrintf(PETSC_COMM_WORLD, "v,5,%f,%f,%f\n",
                      tmp(0), tmp(1), tmp(2));

          PetscPrintf(PETSC_COMM_WORLD, "v,6,%f,%f,%f\n",
                      P0.getv()(0), P0.getv()(1), P0.getv()(2));
          PetscPrintf(PETSC_COMM_WORLD, "v,7,%f,%f,%f\n\n",
                      P1(0), P1(1), P1(2));
        }
          

        // If we are inside the ib, find the location where the particle
        // path intersected the surface.
        if (dir >= GeomX::Tol(0.0))
        {
          double dist;
          Vector3d PN;
          if (! raytrace(P0, P1, dist, user->rtreeSurfs, result, user, pId))
          {
            user->globalError |= SURF_INTERS_PART_TRAJ;
            PetscPrintf(PETSC_COMM_WORLD,
                        "We've got a bit of a problem "
                        "in surface intersection in ParticleTrajectory!\n");
          }
          else
          {
            if (wroDebugPID == -1 || wroDebugPID == pId)
            {
              PetscPrintf(PETSC_COMM_WORLD, "\nParticleTrajectory\n");
              PetscPrintf(PETSC_COMM_WORLD, "v,8,%f,%f,%f\n\n",
                                                           P1(0), P1(1), P1(2));
            }


            result = user->rtreeGrid->getClosestEntity(P1, &dist, searchRad);
            if (result.isValid() && dist <= GeomX::Tol(0.0))
            {
              Id(result).ijk(loc1.i, loc1.j, loc1.k);
              updatePartStatus(pId, ON_IB_SURFACE, loc1, P1, P);
              P->geomLookupStatus = GEOM_SUCCESSFUL;
            }
            else
            {
              user->globalError |= IB_GRID_LOOKUP_ERROR;
              PetscPrintf(PETSC_COMM_WORLD,
                          "Error in grid lookup in ib routine.\n");
            }
          }
        }
      }      
    }
  }


  P->Vp.set(Vp1);

  P->loc = loc1;
  P->coor.x = P1(0);
  P->coor.y = P1(1);
  P->coor.z = P1(2);
  // P->ztraj.push_back(P->coor);
    
  if (user->useDebugDrag)
  {
    user->debugDrag->add(user, pId, time, Vf0.get(), Vp0.get(), Vpf0, P->DDiaPt,
                         user->DDenFd, P->DDenPt, DMasPt, DAccGrav(2),
                         rey, P->RespTime, Cd, DFdrag1, DFdrag2, DFdrag3,
                         DFdrag, DFtot, DAcc, DAcc * (T*ptc_dt) / U + Vpf0,
                         Vt, Vf1.get(), Vp1, Vpf1, P0, loc0, P1, loc1, P->partStat);
  }


  // pId time vp pc vp pc
  // if (pId >= 0 && pId < 100)
  // if (user->RKtest && (pId == 0 || pId == 45 || pId == 91))
  if (pId == 0 || pId == 1 || pId == 9 || pId == 45 || pId == 91 || pId == 98)
  {
    // double err = (Vp1-RKresult.getv()).mag() + (P1-RKresult.getx()).mag();
    // fprintf(user->rkdebug, "%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
    //         pId, time, Vp1(0), Vp1(1), Vp1(2), Vp1.mag(), P1(0), P1(1), P1(2),
    //         RKresult.getv()(0), RKresult.getv()(1), RKresult.getv()(2), RKresult.getv().mag(),
    //         RKresult.getx()(0), RKresult.getx()(1), RKresult.getx()(2));

    fprintf(user->rkdebug, "%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%d\n",
            pId, time, Vp1(0), Vp1(1), Vp1(2), Vp1.mag(), P1(0), P1(1), P1(2), user->calcMethod);
  }


  // DEBUGGING - is the timestep too big?
  if (P->geomLookupStatus == GEOM_SUCCESSFUL)
  {
    int maxDelta = std::max(abs(loc0.i - P->loc.i), abs(loc0.j - P->loc.j));
        maxDelta = std::max(maxDelta,               abs(loc0.k - P->loc.k));
    if (maxDelta > 1)
    {
      user->globalError |= MOVED_TOO_FAR;
      PetscPrintf(PETSC_COMM_WORLD,
                  "pId: %d  Timestep is too big!  delta: %d\n", pId, maxDelta);
    }
  }

  return 0;
}


void updatePartStatus(int pId, ParticleStatus partStat, Location const & loc,
                      Vector3d const & coor, Particle *P)
{
  static char const * pfs =
    "Particle %d at i,j,k: %d,%d,%d  x,y,z: %f,%f,%f  status %d (%s).\n";

  PetscPrintf(PETSC_COMM_WORLD, pfs, pId, loc.i, loc.j, loc.k,
              coor(0), coor(1), coor(2), partStat, partMess[partStat]);
  P->partStat = partStat;
  if (P->partStat == IN_IB_CELL       ||
      P->partStat == IN_SOLID_CELL    ||
      P->partStat == OUTSIDE_AT_INLET ||
      P->partStat == OUTSIDE_AT_WALL   ||
      P->partStat == ON_IB_SURFACE)
    P->Vf0 = P->Vf1 = P->Vp = Vector3d(0.0);
}


bool raytrace(Vector3d const & P0, Vector3d & P1, double & dist,
              std::unique_ptr<Geom::Rtree> & rtree, Geom::Id & result,
              std::shared_ptr< UserCtx > & user, int pId)
{
  auto dir = P0 - P1;
  dir.normalize();
  Ray3Track ray(P1, dir, user, pId);
  rtree->beginRayQuery(ray);
  result = rtree->query(&dist);
  rtree->endQuery();
  if (result.isValid())
    P1 = ray.getOrigin() + ray.getDirection() * dist;
  return result.isValid();
}


// If the particle is outside the grid, find the intersection point
// with the face along P0->P1.
//
// 2022-01-03
// There is an error by searching from P1 to P0. The ray does always
// intersect with the selected face. Going from the inside out will
// always find a face and the search can continue if the found face
// is a face inside the mesh.

// RAY_TRACE_TYPE: 0 = original code
//                 1 = new code
#define RAY_TRACE_TYPE 1

void raytrace2(Vector3d const & P0, Vector3d & P1,
               int pId, Particle * P, std::shared_ptr< UserCtx > & user,
               Geom::Id & result, Location & loc, double searchRad)
{
  // Do the fast search first by limiting the search distance. This is
  // the fastest way and will find the cell 99% of the time.
  double dist;

  // The value of this variable is the external face id of the cell
  // using the getCellSidedness below.
  int bndyFaceId = -1;

  result = user->rtreeGrid->getClosestEntity(P1, &dist, searchRad);
  if (result.isValid())
  {
    // We've found a cell, we are going to need its index location.
    Id(result).ijk(loc.i, loc.j, loc.k);
    P->geomLookupStatus = GEOM_SUCCESSFUL;

    // getClosesEntity find the closest cell that the location P1 it in.
    // It returns a dist of zero if the point in inside or on the structured
    // mesh outer boundary. However, it can also find cells even if the
    // point is outside the mesh.  First check, if the distance is > zero,
    // it MUST be outside the mesh.  In this case, we need to go the part
    // of this method that calculates the exit point.
    if (dist <= GeomX::Tol(0.0))
    {
      // For this next part, there are two types of cells: cell that have
      // a face on the outside boundary of the structure mesh and the
      // other cells which are internal to the mesh.  We do a quick check
      // to see if the cell is internal and, if it is, we're done and
      // can return. If the cell is on the boundary, the routine returns
      // the id of the boundary face.
      if (SpatialIndexGrid::getCellSidedness(loc.i, loc.j, loc.k, bndyFaceId) ==
                                                    SpatialIndexGrid::INTERNAL)
        return;

      // We now know that the cell is on the structured mesh boundary.
      // We need to know if it is inside, outside or on a face of the cell.
      // The next check will only give a definitive answer if it is outside.
      // If it is outside, we need to go the part of this method that
      // calculates the exit point otherwise we need to check to see if
      // the point is internal to the cell are on the external boundary.
      Cmpnts p, coef;
      p.x = P1(0);
      p.y = P1(1);
      p.z = P1(2);
      if (ISPointInCell(p, user.get(), loc, &coef, P, pId, false))
      {
        // Check to see if the point is internal to the cell or on the
        // external boundary. If it is internal, we're done.
        std::vector<Vector3d> fcoor(4);
        SpatialIndexGrid::getFaceCoords(bndyFaceId, loc.i, loc.j, loc.k,
                                        user->coor, fcoor);
        if (sqrt(Geom::facePointDistanceSquared(fcoor, P1)) != GeomX::Tol(0.0))
          return;
      }
    }
  }
  // If the above search doesn't find a cell, it must be outside so we
  // don't need to do anything else since the following code is
  // designed to find the intersection point. This comment is here
  // because which did a getClosestEntity with a maximum search radius
  // the could be very slow and wasn't needed anyway.


  // Need to do another ray query to find exit point and identify which side
  // of the grid the particle exited.
#if RAY_TRACE_TYPE == 0
  auto dir = P0 - P1;
  dir.normalize();
  Ray3Track ray(P1, dir);
#else
  auto dir = P1 - P0;
  dir.normalize();
  Ray3Track ray(P0, dir, user, pId);
  result = user->rtreeGrid->getClosestEntity(P0, &dist, searchRad);
#endif
  bool OK = user->spatialIndexGrid->intersectsRay(result, ray, &dist);
  if (! OK)
  {
    user->globalError |= RAY_QUERY_RAYTRACE2;
    PetscPrintf(PETSC_COMM_WORLD,"PROBLEM: Ray query in raytrace2 failed.\n");
    updatePartStatus(pId, OUTSIDE_UNKNOWN, loc, P1, P);
    P->geomLookupStatus = GEOM_FAILED;
  }
  else
  {
    // Make sure the exit path goes throught the correct face.
    auto fId = ray.getFaceId();

    if (bndyFaceId != -1 && fId != -1 && bndyFaceId != fId)
      return;

    P1 = ray.getOrigin() + ray.getDirection() * dist;
    
    if (fId == -1)
    {
      updatePartStatus(pId, OUTSIDE_UNKNOWN, loc, P1, P);
    }
    else
    {
      // face id   desc
      //    0      max i
      //    1      min i
      //    2      max j
      //    3      min j
      //    4      max k
      //    5      min k
      if (fId == 4)
      {
        P->partStat = OUTSIDE_AT_OUTLET;
      }
      else if (fId == 5)
      {
        P->partStat = OUTSIDE_AT_INLET;
      }
      else if (fId == -1)
      {
        user->globalError |= FACE_ID_RAYTRACE2;
        PetscPrintf(PETSC_COMM_WORLD,
                                 "PROBLEM retreiving face Id in raytrace2.\n");
        P->partStat = OUTSIDE_UNKNOWN;
      }
      else
      {
        P->partStat = OUTSIDE_AT_WALL;
      }

      updatePartStatus(pId, P->partStat, loc, P1, P);
    }
  }
}


PetscErrorCode VelocityReadin(std::shared_ptr< UserCtx > & user, PetscInt ti, PetscInt bi)
{
  PetscViewer viewer;

  char filen[90];

  VecCopy(user->Ucat, user->Ucat_o);

  //sprintf(filen, "ufield%06d_%1d.dat", ti, bi);
  sprintf(filen, "ufield%06d_0.dat", ti);

  // check that the file exists
  auto fd = fopen(filen, "r");
  if (!fd)
  {
    PetscPrintf(PETSC_COMM_WORLD, "Can't open file %s.\n", filen);
    return 1;
  }
  fclose(fd);


  int ierr =
       PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  if (ierr)
  {
    PetscPrintf(PETSC_COMM_WORLD, "Error opening the file %s\n", filen);
    return ierr;
  }

  VecLoadIntoVector(viewer, user->Ucat_c);
  PetscViewerDestroy(viewer);

  PetscInt i, j, k;

  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs;
  PetscInt  	ys = info.ys;
  PetscInt	zs = info.zs;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  Cmpnts ***ucat, ***ucat_c;
  DAVecGetArray(user->fda, user->Ucat, &ucat);
  DAVecGetArray(user->fda, user->Ucat_c, &ucat_c);
  for (k=zs; k<mz-1; k++)
  for (j=ys; j<my-1; j++)
  for (i=xs; i<mx-1; i++) {
    if (i==0 || j==0 || k==0 || i==mx-2 || j==my-2 || k==mz-2)
    {
      ucat[k][j][i].x = 0.;
      ucat[k][j][i].y = 0.;
      ucat[k][j][i].z = 0.;
    }
    else
      ucat[k][j][i] = CENT(ucat_c, k, j, i);
  }
  DAVecRestoreArray(user->fda, user->Ucat, &ucat);
  DAVecRestoreArray(user->fda, user->Ucat_c, &ucat_c);

  return 0;
}


template<typename T>
void vtkOutputLoop(FILE *fd, int cnt, int nPts, double time,
                   char const * fmt, std::vector<Particle> & Ptc,
                   std::function<T(std::vector<Particle> &, int, double)> get)
{
  int cnt2 = 0;
  for (int pn=0; pn<nPts; pn+=1)
  {
    if (Ptc[pn].stime <= GeomX::Tol(time))
    {
      ++cnt2;
      if (cnt2 > cnt)
        break;

      fprintf(fd, fmt, get(Ptc, pn, time));
    }
  }
}


template<typename T>
void vtkOutputScalar(FILE *fd, int cnt, int nPts, double time, char const * label,
                     char const * fmt, std::vector<Particle> & Ptc,
                     std::function<T(std::vector<Particle> &, int, double)> get)
{
  fprintf(fd, "POINT_DATA %d\n", cnt);
  fprintf(fd, "SCALARS %s float 1\n", label);

  vtkOutputLoop(fd, cnt, nPts, time, fmt, Ptc, get);
}


double vtkGetU       (std::vector<Particle> & Ptc, int idx, double time) { return Ptc[idx].Vp.getv()(0); }
double vtkGetV       (std::vector<Particle> & Ptc, int idx, double time) { return Ptc[idx].Vp.getv()(1); }
double vtkGetW       (std::vector<Particle> & Ptc, int idx, double time) { return Ptc[idx].Vp.getv()(2); }
double vtkGetSTIME   (std::vector<Particle> & Ptc, int idx, double time) { return Ptc[idx].stime; }
double vtkGetETIME   (std::vector<Particle> & Ptc, int idx, double time) { return Ptc[idx].etime; }
double vtkGetDIA     (std::vector<Particle> & Ptc, int idx, double time) { return Ptc[idx].DDiaPt; }
double vtkGetDEN     (std::vector<Particle> & Ptc, int idx, double time) { return Ptc[idx].DDenPt; }
int    vtkGetSTATUS  (std::vector<Particle> & Ptc, int idx, double time) { return Ptc[idx].partStat; }
double vtkGetTIME    (std::vector<Particle> & Ptc, int idx, double time) { return time; }
int    vtkGetPID     (std::vector<Particle> & Ptc, int idx, double time) { return idx; }
double vtkGetREY     (std::vector<Particle> & Ptc, int idx, double time) { return Ptc[idx].Rey; }
double vtkGetRESPTIME(std::vector<Particle> & Ptc, int idx, double time) { return Ptc[idx].RespTime; }
int    vtkGetBID     (std::vector<Particle> & Ptc, int idx, double time) { return floor(Ptc[idx].stime / 5.0) + 1; }


PetscErrorCode ParticleTrajectoryvtk(std::shared_ptr< UserCtx > & user,
                                     std::vector<Particle> & Ptc,
                                     PetscInt timestep, double time)
{
  FILE *fd;
  char filen[80];
  sprintf(filen, "Particle.%06d.vtk", timestep, 0);
  fd = fopen(filen, "w");

  fprintf(fd, "# vtk DataFile Version 2.0\n");
  fprintf(fd, "Particle positions at time %f\n", time);
  fprintf(fd, "ASCII\n");
  fprintf(fd, "DATASET UNSTRUCTURED_GRID\n");


  // Count number of points still inside the grid.
  int cnt = 0;
  for (int pn=0; pn<user->ptc_num; ++pn)
  {
    if (Ptc[pn].stime <= GeomX::Tol(time))
      ++cnt;
  }
  fprintf(fd, "POINTS %d float\n", cnt);


  int cnt2 = 0;
  for (int pn=0; pn<user->ptc_num; pn+=1)
  {
    if (Ptc[pn].stime <= GeomX::Tol(time))
    {
      ++cnt2;
      if (cnt2 > cnt)
        break;

      fprintf(fd, "%.6e %.6e %.6e\n", Ptc[pn].coor.x, Ptc[pn].coor.y, Ptc[pn].coor.z);
    }
  }


  cnt2 = 0;
  fprintf(fd, "CELLS %d %d\n", cnt, 2*cnt);
  for (int pn=0; pn<user->ptc_num; pn+=1)
  {
    if (Ptc[pn].stime <= GeomX::Tol(time))
    {
      ++cnt2;
      if (cnt2 > cnt)
        break;

      fprintf(fd, "1 %d\n", cnt2);
    }
  }

  cnt2 = 0;
  fprintf(fd, "CELLS_TYPES %d\n", cnt);
  for (int pn=0; pn<user->ptc_num; pn+=1)
  {
    if (Ptc[pn].stime <= GeomX::Tol(time))
    {
      ++cnt2;
      if (cnt2 > cnt)
        break;

      fprintf(fd, "1\n");
    }
  }

  vtkOutputScalar<double>(fd, cnt, user->ptc_num, time, "U", "%.5e\n",        Ptc, vtkGetU);
  vtkOutputScalar<double>(fd, cnt, user->ptc_num, time, "V", "%.5e\n",        Ptc, vtkGetV);
  vtkOutputScalar<double>(fd, cnt, user->ptc_num, time, "W", "%.5e\n",        Ptc, vtkGetW);
  vtkOutputScalar<double>(fd, cnt, user->ptc_num, time, "STIME", "%.5e\n",    Ptc, vtkGetSTIME);
  vtkOutputScalar<double>(fd, cnt, user->ptc_num, time, "ETIME", "%.5e\n",    Ptc, vtkGetETIME);
  vtkOutputScalar<double>(fd, cnt, user->ptc_num, time, "DIA", "%.5e\n",      Ptc, vtkGetDIA);
  vtkOutputScalar<double>(fd, cnt, user->ptc_num, time, "DEN", "%.5e\n",      Ptc, vtkGetDEN);
  vtkOutputScalar<int>   (fd, cnt, user->ptc_num, time, "STATUS", "%d\n",     Ptc, vtkGetSTATUS);
  vtkOutputScalar<double>(fd, cnt, user->ptc_num, time, "TIME", "%.5e\n",     Ptc, vtkGetTIME);
  vtkOutputScalar<int>   (fd, cnt, user->ptc_num, time, "PID", "%d\n",        Ptc, vtkGetPID);
  vtkOutputScalar<double>(fd, cnt, user->ptc_num, time, "REY", "%.5e\n",      Ptc, vtkGetREY);
  vtkOutputScalar<double>(fd, cnt, user->ptc_num, time, "RESPTIME", "%.5e\n", Ptc, vtkGetRESPTIME);
  if (user->useBreath == 1)
    vtkOutputScalar<int>( fd, cnt, user->ptc_num, time, "BID", "%d\n",        Ptc, vtkGetBID);

  fclose(fd);
  return 0;
}


PetscErrorCode ParticleTrajectoryIO2(std::shared_ptr< UserCtx > & user,
                                     std::vector<Particle> & Ptc,
                                     PetscInt timestep, double time)
{
  auto elapTime = user->globTimer.getElapsedCPUTimeStr();
  PetscPrintf(PETSC_COMM_WORLD, "Writing particle file %d for time %f.  "
              "Elapsed time: %s\n",
              timestep, time, elapTime.c_str());


  if (user->vtk)
    return ParticleTrajectoryvtk(user, Ptc, timestep, time);


  FILE *fd;
  char filen[80];  
  // sprintf(filen, "Particle%06d_%1d.dat", timestep, 0);
  sprintf(filen, "Particle.%06d.dat", timestep);
  fd = fopen(filen, "w");

  // PetscFPrintf(PETSC_COMM_WORLD, fd, "Variables = X, Y, Z, Out\n");
  if (user->useBreath == 1)
    PetscFPrintf(PETSC_COMM_WORLD, fd,
                 "Variables =\tX,\tY,\tZ,\tU,\tV,\tW,"
                 "\tSTIME,\tETIME,\tDIA,\tDEN,\tOUT,\tTIME,\tPID,\tBID,\tREY,\tRESPTIME\n");
  else
    PetscFPrintf(PETSC_COMM_WORLD, fd,
                 "Variables =\tX,\tY,\tZ,\tU,\tV,\tW,"
                 "\tSTIME,\tETIME,\tDIA,\tDEN,\tOUT,\tTIME,\tPID,\tREY,\tRESPTIME\n");

  // Count number of points still inside the grid.
  int cnt = 0;
  for (int pn=0; pn<user->ptc_num; ++pn)
  {
    if (Ptc[pn].stime <= GeomX::Tol(time))
      ++cnt;
  }

  PetscInt pn;
  if (cnt == 0)
  {
    PetscFPrintf(PETSC_COMM_WORLD, fd,
                 "ZONE I = %i,\tJ = %i,\tDATAPACKING = POINT\n", 1, 1);
    PetscFPrintf(PETSC_COMM_WORLD, fd, "1000.0\t1000.0\t1000.0\t1\n");
  }
  else
  {
    PetscFPrintf(PETSC_COMM_WORLD, fd,
                 "ZONE I = %i,\tJ = %i,\tDATAPACKING = POINT\n", 1, cnt);
    int cnt2 = 0;
    for (pn=0; pn<user->ptc_num; pn+=1)
    {
      if (Ptc[pn].stime <= GeomX::Tol(time))
      {
        ++cnt2;
        if (cnt2 > cnt)
          break;

        if (user->useBreath == 1)
        {
          int bId = floor(Ptc[pn].stime / 5.0) + 1;
          PetscFPrintf(PETSC_COMM_WORLD, fd,
               "%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%d\t%le\t%d\t%d\t%le\t%le\n",
                   Ptc[pn].coor.x, Ptc[pn].coor.y, Ptc[pn].coor.z,
                   Ptc[pn].Vp.getv()(0), Ptc[pn].Vp.getv()(1), Ptc[pn].Vp.getv()(2), Ptc[pn].stime,
                   Ptc[pn].etime, Ptc[pn].DDiaPt, Ptc[pn].DDenPt,
                   Ptc[pn].partStat, time, pn, bId, Ptc[pn].Rey, Ptc[pn].RespTime);
        }
        else
          PetscFPrintf(PETSC_COMM_WORLD, fd,
               "%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%d\t%le\t%d\t%le\t%le\n",
                   Ptc[pn].coor.x, Ptc[pn].coor.y, Ptc[pn].coor.z,
                   Ptc[pn].Vp.getv()(0), Ptc[pn].Vp.getv()(1), Ptc[pn].Vp.getv()(2), Ptc[pn].stime,
                   Ptc[pn].etime, Ptc[pn].DDiaPt, Ptc[pn].DDenPt,
                   Ptc[pn].partStat, time, pn, Ptc[pn].Rey, Ptc[pn].RespTime);
      }
    }
  }

  fclose(fd);
  return 0;
}



// Write particle trajectory tracks in new format
// timestep is the index to be used for the filename.
// time is the time of ti.
PetscErrorCode ParticleTrajectoryIO3(std::shared_ptr< UserCtx > & user,
                                     std::vector<Particle> & Ptc,
                                     PetscInt timestep, double time)
{
  // This is really bad coding!  Don't do it!
  static int currTimeStep = -1;

  FILE * fd;
  char filename[128];
  sprintf(filename, "Part.Traj.%06d.dat", timestep);
  
  if (currTimeStep == timestep)
  {
    fd = fopen(filename, "a");
  }
  else
  {
    fd = fopen(filename ,"w");
    currTimeStep = timestep;
    PetscFPrintf(PETSC_COMM_WORLD, fd,
              "time\tpId\tcoor.x\tcoor.y\tcoor.z\tvel.x\tvel.y\tvel.z\tstime\t"
              "etime\tdia\tden\tpartStat\n");
  }

  for (int pn=0; pn<user->ptc_num; pn++)
  {
    if (Ptc[pn].stime <= GeomX::Tol(time))
      PetscFPrintf(PETSC_COMM_WORLD, fd,
             "%le\t%d\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%d\n",
                   time, pn, Ptc[pn].coor.x, Ptc[pn].coor.y, Ptc[pn].coor.z,
                   Ptc[pn].Vp.getv()(0), Ptc[pn].Vp.getv()(1), Ptc[pn].Vp.getv()(2), Ptc[pn].stime,
                   Ptc[pn].etime, Ptc[pn].DDiaPt, Ptc[pn].DDenPt, Ptc[pn].partStat);
  }

  fclose(fd);
}



PetscErrorCode BoundingSphere(std::shared_ptr< UserCtx > & user)
{
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs;
  PetscInt  	ys = info.ys;
  PetscInt	zs = info.zs;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt i, j, k;

  Vec Coor;
  PetscScalar ***rad;
  Cmpnts ***coor, ***cent;
  DAGetCoordinates(user->da, &Coor);
  DAVecGetArray(user->fda, Coor, &coor);
  DAVecGetArray(user->fda, user->Cent, &cent);
  DAVecGetArray(user->fda, user->Rad, &rad);

  Cmpnts v;
  VecSet(user->Rad, 0.0);

  for (k=zs; k<mz-2; k++) {
    for (j=ys; j<my-2; j++) {
      for (i=xs; i<mx-2; i++) {
	VecAB(v, cent[k][j][i], coor[k][j][i]);
	rad[k][j][i] = PetscMax(rad[k][j][i], Length(v));
	VecAB(v, cent[k][j][i], coor[k][j][i+1]);
	rad[k][j][i] = PetscMax(rad[k][j][i], Length(v));
	VecAB(v, cent[k][j][i], coor[k][j+1][i]);
	rad[k][j][i] = PetscMax(rad[k][j][i], Length(v));
	VecAB(v, cent[k][j][i], coor[k+1][j][i]);
	rad[k][j][i] = PetscMax(rad[k][j][i], Length(v));
	VecAB(v, cent[k][j][i], coor[k][j+1][i+1]);
	rad[k][j][i] = PetscMax(rad[k][j][i], Length(v));
	VecAB(v, cent[k][j][i], coor[k+1][j][i+1]);
	rad[k][j][i] = PetscMax(rad[k][j][i], Length(v));
	VecAB(v, cent[k][j][i], coor[k+1][j+1][i]);
	rad[k][j][i] = PetscMax(rad[k][j][i], Length(v));
	VecAB(v, cent[k][j][i], coor[k+1][j+1][i+1]);
	rad[k][j][i] = PetscMax(rad[k][j][i], Length(v));

	rad[k][j][i] += 0.001;
      }
    }
  }

  DAVecRestoreArray(user->fda, Coor, &coor);
  DAVecRestoreArray(user->fda, user->Cent, &cent);
  DAVecRestoreArray(user->fda, user->Rad, &rad);

  return 0;
}

typedef struct {
  PetscReal xmin, xmax, ymin, ymax, zmin, zmax;

  PetscInt IM, JM, KM;

  PetscReal dx, dy, dz;
} CoorBound;


CoorBound gridbound;

void PointLocation(Location & L, Cmpnts & P)
{
  L.i = (P.x - gridbound.xmin) / gridbound.dx;
  L.j = (P.y - gridbound.ymin) / gridbound.dy;
  L.k = (P.z - gridbound.zmin) / gridbound.dz;
}

#define CellLoc(i, j, k) \
  k * gridbound.JM*gridbound.IM + j * gridbound.IM + i

PetscErrorCode UniformGrid(std::shared_ptr< UserCtx > & user)
{
  PetscPrintf(PETSC_COMM_WORLD, "DEBUG entering UniformGrid\n");
  Util::Timer timer(Util::Timer::Wall);

  PetscInt i, j, k;

  PetscInt mi, mj, mk;

  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;

  gridbound.IM = 40; gridbound.JM = 40; gridbound.KM = 40;

  gridbound.xmin = user->coor[0][0][0].x; gridbound.xmax = gridbound.xmin;
  gridbound.ymin = user->coor[0][0][0].y; gridbound.ymax = gridbound.ymin;
  gridbound.zmin = user->coor[0][0][0].z; gridbound.zmax = gridbound.zmin;

  PetscMalloc(gridbound.IM*gridbound.JM*gridbound.KM*sizeof(List),
	      &user->CellList);

  PetscReal tmp;

  for (k=0; k<gridbound.KM; k++) {
    for (j=0; j<gridbound.JM; j++) {
      for (i=0; i<gridbound.IM; i++) {
	initlist(&user->CellList[CellLoc(i, j, k)]);
      }
    }
  }

  for (k=zs; k<ze-1; k++) {
    for (j=ys; j<ye-1; j++) {
      for (i=xs; i<xe-1; i++) {
	if ((tmp = user->coor[k][j][i].x) < gridbound.xmin) {
	  gridbound.xmin = tmp;
	}
	else if (tmp > gridbound.xmax) {
	  gridbound.xmax = tmp;
	}
	if ((tmp = user->coor[k][j][i].y) < gridbound.ymin) {
	  gridbound.ymin = tmp;
	}
	else if (tmp > gridbound.ymax) {
	  gridbound.ymax = tmp;
	}
	if ((tmp = user->coor[k][j][i].z) < gridbound.zmin) {
	  gridbound.zmin = tmp;
	}
	else if (tmp > gridbound.zmax) {
	  gridbound.zmax = tmp;
	}
      }
    }
  }

  gridbound.dx = (gridbound.xmax - gridbound.xmin) / (gridbound.IM -1);
  gridbound.dy = (gridbound.ymax - gridbound.ymin) / (gridbound.JM -1);
  gridbound.dz = (gridbound.zmax - gridbound.zmin) / (gridbound.KM -1);

  PetscInt iradial, jradial, kradial;
  PetscInt mis, mie, mjs, mje, mks, mke;

  Location sloc, cell;

  for (k=zs; k<ze-2; k++) {
    // PetscPrintf(PETSC_COMM_WORLD, "Insert node at plane %i\n", k);
    for (j=ys; j<ye-2; j++) {
      for (i=xs; i<xe-2; i++) {
	PointLocation(sloc, user->cent[k][j][i]);
	iradial = user->rad[k][j][i] / gridbound.dx;
	jradial = user->rad[k][j][i] / gridbound.dy;
	kradial = user->rad[k][j][i] / gridbound.dz;

	mis = PetscMax(0, sloc.i - iradial);
	mie = PetscMin(gridbound.IM, sloc.i + iradial + 1);
	mjs = PetscMax(0, sloc.j - jradial);
	mje = PetscMin(gridbound.JM, sloc.j + jradial + 1);
	mks = PetscMax(0, sloc.k - kradial);
	mke = PetscMin(gridbound.KM, sloc.k + kradial + 1);

	cell.i = i; cell.j = j; cell.k = k;
	for (mk = mks; mk < mke; mk++) {
	  for (mj = mjs; mj < mje; mj++) {
	    for (mi = mis; mi < mie; mi++) {
	      insertnode(&user->CellList[CellLoc(mi, mj, mk)], cell);
	    }
	  }
	}
      }
    }
  }

  std::stringstream mess;
  mess << timer;
  PetscPrintf(PETSC_COMM_WORLD, "%s\n", mess.str().c_str());
  PetscPrintf(PETSC_COMM_WORLD, "DEBUG exiting UniformGrid\n");

  return 0;
}


PetscErrorCode ParticleLocation1(std::shared_ptr< UserCtx > & user, Particle * P, int pId)
{
  Cmpnts p = P->coor;

  Location searchingcell;

  PointLocation(searchingcell, p);
  PetscTruth NotDecided = PETSC_TRUE;

  node *current;
  PetscTruth ISInside;
  while (NotDecided) {
    current = user->CellList[CellLoc(searchingcell.i, searchingcell.j, searchingcell.k)].head;

    while (current) {
/*       PetscPrintf(PETSC_COMM_WORLD, "Searching in cell %i %i %i\n", current->Node.i, current->Node.j, current->Node.k); */
      if ( (ISInside = ISPointInCell(p, user.get(), current->Node, &P->coef, P, pId)) )
      {
	NotDecided = PETSC_FALSE;
	P->loc = current->Node;
	current = PETSC_NULL;
	break;
      }
      else {
	current = current->next;
      }
    }
    if (!current) break;
  }

  if (NotDecided) {
    current = current;
  }
  return (!NotDecided);
}


void PrintDebugInfo01(bool header, Particle *P, int pId, int ti,
                      Vector3d const & DAcc, Vector3d const & Vf,
                      Vector3d const & Vpf1, Vector3d const & Vp1,
                      Vector3d const & P1)
{
  if (wroDebugPID != -1 && wroDebugPID != pId)
    return;

  FILE * fd;
  char filename[128];
  char const * mode = (header ? "w" : "a");

  sprintf(filename, "debug_acc_%06d.dat", pId);
  fd = fopen(filename, mode);
  if (header)
    fprintf(fd, "Time Index\tDAcc-x\tDAcc-y\tDAcc-z\tDAcc-mag\n");
  else
    fprintf(fd, "%d\t%e\t%e\t%e\t%e\n", ti, DAcc(0), DAcc(1), DAcc(2),
            DAcc.mag());
  fclose(fd);

    
  sprintf(filename, "debug_Vf_%06d.dat", pId);
  fd = fopen(filename, mode);
  if (header)
    fprintf(fd, "Time Index\tVf-x\tVf-y\tVf-z\tVf-mag\n");
  else
    fprintf(fd, "%d\t%e\t%e\t%e\t%e\n", ti, Vf(0), Vf(1), Vf(2), Vf.mag());
  fclose(fd);

    
  sprintf(filename, "debug_Vpf1_%06d.dat", pId);
  fd = fopen(filename, mode);
  if (header)
    fprintf(fd, "Time Index\tVpf1-x\tVpf1-y\tVpf1-z\tVpf1-mag\n");
  else
    fprintf(fd, "%d\t%e\t%e\t%e\t%e\n", ti,
                                         Vpf1(0), Vpf1(1), Vpf1(2),Vpf1.mag());
  fclose(fd);

    
  sprintf(filename, "debug_Vp1_%06d.dat", pId);
  fd = fopen(filename, mode);
  if (header)
    fprintf(fd, "Time Index\tVp1-x\tVp1-y\tVp1-z\tVp1-mag\n");
  else
    fprintf(fd, "%d\t%e\t%e\t%e\t%e\n", ti, Vp1(0), Vp1(1),Vp1(2),Vp1.mag());
  fclose(fd);

    
  sprintf(filename, "debug_P1_%06d.dat", pId);
  fd = fopen(filename, mode);
  if (header)
    fprintf(fd, "Time Index\tP1-x\tP1-y\tP1-z\tP1-mag\n");
  else
    fprintf(fd, "%d\t%e\t%e\t%e\t%e\n", ti, P1(0), P1(1), P1(2),P1.mag());
  fclose(fd);
}


void PrintDebugInfo02(int pId, Vector3d const & P0, double ptc_dt,
                      double DVolPt, double DMasPt, Vector3d const & DAccGrav,
                      Vector3d const & Vf, Vector3d const & Vp0,
                      Vector3d const & Vpf0, double rey, double Cd,
                      Vector3d const & DFdrag, Vector3d const & DFtot,
                      Vector3d const & DAcc, Vector3d const & DVel,
                      Vector3d const & Vpf1, double const & Vt,
                      Vector3d const & Vp1, Vector3d const & Pdel,
                      Vector3d const & P1)
{
  if (wroDebugPID != -1 && wroDebugPID != pId)
    return;

  PetscPrintf(PETSC_COMM_WORLD,
              "\nP0: %f,%f,%f\n", P0(0), P0(1), P0(2));
  PetscPrintf(PETSC_COMM_WORLD, "ptc_dt: %f\n", ptc_dt);
  PetscPrintf(PETSC_COMM_WORLD, "DVolPt: %e\n", DVolPt);
  PetscPrintf(PETSC_COMM_WORLD, "DMasPt: %e\n", DMasPt);
  PetscPrintf(PETSC_COMM_WORLD,
              "DAccGrav: %e,%e,%e\n", DAccGrav(0), DAccGrav(1), DAccGrav(2));
  PetscPrintf(PETSC_COMM_WORLD,
              "Vp0: %f,%f,%f  mag: %f\n", Vp0(0), Vp0(1), Vp0(2), Vp0.mag());
  PetscPrintf(PETSC_COMM_WORLD,
              "Vf: %f,%f,%f  mag: %f\n", Vf(0), Vf(1), Vf(2), Vf.mag());
  PetscPrintf(PETSC_COMM_WORLD,
              "Vpf0: %f,%f,%f  mag: %f\n",
                                        Vpf0(0), Vpf0(1), Vpf0(2), Vpf0.mag());
  PetscPrintf(PETSC_COMM_WORLD, "Vratio: %f\n", Vp0.mag() / Vf.mag());
  PetscPrintf(PETSC_COMM_WORLD, "rey: %f\n", rey);
  PetscPrintf(PETSC_COMM_WORLD, "Cd: %f\n", Cd);
  PetscPrintf(PETSC_COMM_WORLD,
              "DFdrag: %e,%e,%e\n", DFdrag(0), DFdrag(1), DFdrag(2));
  PetscPrintf(PETSC_COMM_WORLD,
              "DFtot: %e,%e,%e\n", DFtot(0), DFtot(1), DFtot(2));
  PetscPrintf(PETSC_COMM_WORLD, "DAcc: %e,%e,%e  mag: %e\n",
              DAcc(0), DAcc(1), DAcc(2), DAcc.mag());
  PetscPrintf(PETSC_COMM_WORLD, "DVel: %f,%f,%f  mag: %e\n",
              DVel(0), DVel(1), DVel(2), DVel.mag());
  PetscPrintf(PETSC_COMM_WORLD, "Vpf1: %f,%f,%f  mag: %f\n",
                                        Vpf1(0), Vpf1(1), Vpf1(2), Vpf1.mag());
  PetscPrintf(PETSC_COMM_WORLD, "Vt: %f\n", Vt);
  PetscPrintf(PETSC_COMM_WORLD,
              "Vp1: %f,%f,%f  mag: %f\n", Vp1(0), Vp1(1), Vp1(2), Vp1.mag());
  PetscPrintf(PETSC_COMM_WORLD,
              "Vp1/Vf mag: %f\n", Vp1.mag() / Vf.mag());
  PetscPrintf(PETSC_COMM_WORLD,
              "Pdel: %f,%f,%f  mag: %f\n",Pdel(0),Pdel(1),Pdel(2),Pdel.mag());
  PetscPrintf(PETSC_COMM_WORLD,
              "P1: %f,%f,%f\n", P1(0), P1(1), P1(2));
  if (Vpf1.mag() > Vf.mag())
    PetscPrintf(PETSC_COMM_WORLD, "Particle velocity is > fluid velocity\n");
}
