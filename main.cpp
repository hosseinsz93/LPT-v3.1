#include <memory>
#include "variables.h"

#include "DragCalcDebug.h"
#include "GeomRtree.h"
#include "SpatialIndexBB.h"
#include "SpatialIndexGrid.h"
#include "SpatialIndexSurf.h"
#include "SpatialIndexSurfs.h"
#include "Surface.h"
#include "RayFaceIntersect.h"

static char help[] = "Testing programming!";

#undef __FUNCT__
#define __FUNCT__ "main"

PetscInt binary_input = 0;
PetscInt block_number;
PetscInt ptc_tn;
PetscReal ptc_dt;
// Number of particle iterations between the iterations of the results files.
PetscInt ptc_ftn;
// Number of particle interations in each particle iteration.
PetscInt ptc_int;
PetscReal cl = 1.;


void OutputGlobalErrors(std::shared_ptr< UserCtx > & user)
{
  if (user->globalError & MOVED_TOO_FAR)
    PetscPrintf(PETSC_COMM_WORLD, "At least one particle was moved too far. "
                                  "Consider reducing the time step.\n");

  if (user->globalError & GRID_QUERY_RAYTRACE2)
    PetscPrintf(PETSC_COMM_WORLD, "Grid query in raytrace2 failed.\n");

  if (user->globalError & RAY_QUERY_RAYTRACE2)
    PetscPrintf(PETSC_COMM_WORLD, "Ray query in raytrace2 failed.\n");

  if (user->globalError & FACE_ID_RAYTRACE2)
    PetscPrintf(PETSC_COMM_WORLD, "Problem retreiving face Id in raytrace2.\n");

  if (user->globalError & SURF_INTERS_PART_TRAJ)
    PetscPrintf(PETSC_COMM_WORLD, "Problem with surface "
                                  "intersection in ParticleTrajectory.\n");

  if (user->globalError & IB_GRID_LOOKUP_ERROR)
    PetscPrintf(PETSC_COMM_WORLD, "Error in grid lookup in ib routine.\n");

  if (user->globalError & RAY_QUERY_PART_TRAJ)
    PetscPrintf(PETSC_COMM_WORLD, "Ray query error in ParticleTrajectory.\n");

  if (user->globalError & COEF_ERROR)
    PetscPrintf(PETSC_COMM_WORLD, "Error of coef calculation in ISPointInCell.\n");

  if (user->globalError & FORCE_REL_VEC_TO_ZERO)
    PetscPrintf(PETSC_COMM_WORLD, "At least one particle's relative velocity "
                "was set to zero because of a small Reynold's number.\n");
}


int main(int argc, char **argv)
{
  PetscInitialize(&argc, &argv, (char *)0, help);

  // Print current version.
  PetscPrintf(PETSC_COMM_WORLD, "LPT v%s\n", progVersion);

  // std::shared_ptr< UserCtx > user(NULL, NoDeleter<UserCtx>());
  std::shared_ptr< UserCtx > user(NULL, [](UserCtx *) { });
  /* PetscInt ti; */
  PetscTruth flag;
  
  PetscOptionsInsertFile(PETSC_COMM_WORLD,"control.lpt.dat", PETSC_TRUE);
  PetscOptionsGetReal(PETSC_NULL, "-cl", &cl, PETSC_NULL);

  PetscOptionsGetInt(PETSC_NULL, "-binary", &binary_input, PETSC_NULL);
  PetscPrintf(PETSC_COMM_WORLD, "-binary %d\n", binary_input);


  {
    PetscPrintf(PETSC_COMM_WORLD, "DEBUG starting to read grid.dat\n");
    Util::Timer timer(Util::Timer::Wall);

    FILE *fd;
    fd = fopen("grid.dat", "r");
    if (!fd)
    {
      PetscPrintf(PETSC_COMM_WORLD, "Can't open file grid.dat.\n");
      return 1;
    }
  
    if (binary_input)
      fread(&block_number, sizeof(int), 1, fd);
    else
      fscanf(fd, "%i\n", &block_number);
    UserCtx * tmp;
    PetscMalloc(block_number*sizeof(UserCtx), &tmp);
    user.reset(tmp);

    // users variable initiation:
    
    // Need to run ctor for C++ objects with are located in variable
    // user where no ctors where calleds during PetscMalloc.
    // Also need to run dtors before user goes out of scope.
    new(&user->spatialIndexGrid)
      std::shared_ptr<SpatialIndexGrid>(nullptr, [](SpatialIndexGrid *) { });
    new(&user->rtreeGrid) std::unique_ptr<Geom::Rtree>;
    new(&user->surfs)
      std::shared_ptr<Surfaces>(nullptr, [](Surfaces *) { });
    new(&user->spatialIndexSurfs)
      std::shared_ptr<SpatialIndexSurfs>(nullptr, [](SpatialIndexSurfs *) { });
    new(&user->rtreeSurfs) std::unique_ptr<Geom::Rtree>;
    new(&user->debugDrag) std::unique_ptr<DragCalcDebug>;
    new(&user->globTimer) Util::Timer(Util::Timer::Wall);


    user->rkdebug = fopen("rkdebug.dat", "w");
    fprintf(user->rkdebug, "pId\ttime\tovx\tovy\tovz\tovmag\topx\topy\topz\tnvx\tnvy\tnvz\tnvmag\tnpx\tnpy\tnpz\terr1\terr3\berr7\n");



    user->globalError = 0;


    int ierr = Initial(fd, user, timer);
    if (ierr)
      return ierr;
  }

  PetscOptionsGetInt(PETSC_NULL, "-tis", &user->tis, &flag);
  PetscPrintf(PETSC_COMM_WORLD, "-tis %i\n", user->tis);
  if (!flag) {
    PetscPrintf(PETSC_COMM_WORLD, "Need to specify the starting time with -tis TIME!\n");
  }
  PetscOptionsGetInt(PETSC_NULL, "-tie", &user->tie, &flag);
  PetscPrintf(PETSC_COMM_WORLD, "-tie %i\n", user->tie);
  if (!flag) {
    PetscPrintf(PETSC_COMM_WORLD, "Need to specify the ending time with -tie TIME!\n");
  }
  PetscOptionsGetInt(PETSC_NULL, "-ts", &user->ts, &flag);
  PetscPrintf(PETSC_COMM_WORLD, "-ts %i\n", user->ts);
  if (!flag) {
    PetscPrintf(PETSC_COMM_WORLD, "Need to specify the time step increasement with -ts TIME_STEP!\n");
  }

  /*
    This option is no longer needed now that the particle information
    is specified in a file.

  PetscOptionsGetInt(PETSC_NULL, "-ptc_num", &user->ptc_num, &flag);
  PetscPrintf(PETSC_COMM_WORLD, "-ptc_num %i\n", user->ptc_num);
  if (!flag) {
    PetscPrintf(PETSC_COMM_WORLD, "Need to specify the number of particles with -ptc_num NUMBER!\n");
   }
  */

  PetscOptionsGetReal(PETSC_NULL, "-ptc_dt", &ptc_dt, &flag);
  PetscPrintf(PETSC_COMM_WORLD, "-ptc_dt %le\n", ptc_dt);
  if (!flag) {
    PetscPrintf(PETSC_COMM_WORLD, "Need to specify particles delta time -ptc_dt!\n");
  }

  PetscOptionsGetReal(PETSC_NULL, "-dt", &user->dt, PETSC_NULL);
  PetscPrintf(PETSC_COMM_WORLD, "-dt %le\n", user->dt);

  
  // Doing a restart  WRO 2021-03-16
  user->ptc_rstart = -1.0;
  PetscOptionsGetReal(PETSC_NULL, "-ptc_rstart", &user->ptc_rstart, PETSC_NULL);
  user->ptc_rstart = (int)user->ptc_rstart;
  if (user->ptc_rstart > 0.0)
    PetscPrintf(PETSC_COMM_WORLD, "We're doing a restart at time step %d\n",
                                                        (int)user->ptc_rstart);
  
  // New velocity calculation  WRO 2021-01-07
  user->calcMethod = 1;
  PetscOptionsGetInt(PETSC_NULL, "-ptc_calc_method", &user->calcMethod, PETSC_NULL);
  PetscPrintf(PETSC_COMM_WORLD, "-ptc_calc_method %d\n", user->calcMethod);

  // Acceration of gravity in x
  user->DAccGrav(0) = 0.0; // m/s^2
  PetscOptionsGetReal(PETSC_NULL, "-ptc_acc_grav-x", &user->DAccGrav(0), PETSC_NULL);

  // Acceration of gravity in y
  user->DAccGrav(1) = 0.0; // m/s^2
  PetscOptionsGetReal(PETSC_NULL, "-ptc_acc_grav-y", &user->DAccGrav(1), PETSC_NULL);

  // Acceration of gravity in z
  user->DAccGrav(2) = 0.0; // m/s^2
  PetscOptionsGetReal(PETSC_NULL, "-ptc_acc_grav-z", &user->DAccGrav(2), PETSC_NULL);

  // Acceration of gravity - for backward compatibility
  user->DAccGrav(2) = 0.0; // m/s^2
  PetscOptionsGetReal(PETSC_NULL, "-ptc_acc_grav", &user->DAccGrav(2), PETSC_NULL);
  user->DAccGrav(2) *= -1.0;

  // Density of fluid - this value is for air at STP
  user->DDenFd = 1.228; // kg/m^3
  PetscOptionsGetReal(PETSC_NULL, "-ptc_den_fluid", &user->DDenFd, PETSC_NULL);

  // Dynamic Viscosity of fluid - this value is for air at 27C
  user->DVisFd = 1.788E-05; // kg/(m s)
  PetscOptionsGetReal(PETSC_NULL, "-ptc_dyn_vis_fluid", &user->DVisFd,PETSC_NULL);

  // Default particle properties
  user->DDiaPt = 30e-6; // m
  PetscOptionsGetReal(PETSC_NULL, "-ptc_part_dia", &user->DDiaPt, PETSC_NULL);

  user->DRoPt = 2500; // kg/m^3
  PetscOptionsGetReal(PETSC_NULL, "-ptc_part_den", &user->DRoPt, PETSC_NULL);

  // For non-dimensional runs, this is length scaling.
  user->lengScale = 1.0; // m
  PetscOptionsGetReal(PETSC_NULL, "-ptc_leng_scale", &user->lengScale, PETSC_NULL);

  // For non-dimensional runs, this is velocity scaling.
  user->velScale = 1.0; // m/s
  PetscOptionsGetReal(PETSC_NULL, "-ptc_vel_scale", &user->velScale, PETSC_NULL);

  // Include ib surfaces as barriers for particle motion.
  user->useIBs = 0;
  PetscOptionsGetReal(PETSC_NULL, "-ptc_use_ibs", &user->useIBs, PETSC_NULL);

  // Method used to initialize particles to be tracked.
  // 0: Don't read ParticleInitial.dat
  // 1: Read from ParticleInitial.dat (default)
  user->initPart = 1;
  PetscOptionsGetReal(PETSC_NULL, "-ptc_init_part", &user->initPart,PETSC_NULL);

  // Use methods to generate/modify particle infomation for Ali's covit
  // breath simulations.  Current, particles are not generated.
  // 1: Use breath routines
  user->useBreath = 0;
  PetscOptionsGetReal(PETSC_NULL, "-ptc_breath", &user->useBreath, PETSC_NULL);

  // If this option is set, the particles' velocity is set to the fluid
  // velocity its the location plus the terminal velocity due to gravity.
  user->setPartVel = 0;
  PetscOptionsGetReal(PETSC_NULL, "-ptc_set_part_vel", &user->setPartVel, PETSC_NULL);

  // Debug drag calculation.
  user->useDebugDrag = 0;
  PetscOptionsGetReal(PETSC_NULL, "-ptc_debug_drag", &user->useDebugDrag, PETSC_NULL);
  if (user->useDebugDrag)
    user->debugDrag.reset(new DragCalcDebug);

  // Debug particle trajectory
  user->useDebugTraj = 0;
  PetscOptionsGetReal(PETSC_NULL, "-ptc_debug_traj", &user->useDebugTraj, PETSC_NULL);

  // Result file output frequency in CFD timestep "units". The default is user->ts.
  user->output_ts = 0;
  PetscOptionsGetInt(PETSC_NULL, "-ptc_output_ts", &user->output_ts, PETSC_NULL);

  // Run the Runge-Kutta integration routines.
  user->RKtest = 0;
  PetscOptionsGetInt(PETSC_NULL, "-RKtest", &user->RKtest, PETSC_NULL);

  // If set, output .vtk files instead of tecplot files.
  user->vtk = 0;
  PetscOptionsGetInt(PETSC_NULL, "-vtk", &user->vtk, PETSC_NULL);

  // ptc_dt must be less than dt and ptc_dt must be an integral multiple of dt.
  // -dt, &user->dt; -ptc_dt, &ptc_dt

  double ratio = user->dt / ptc_dt;
  double fracpart, intpart;
  fracpart = abs(modf(ratio, &intpart));
  if (ptc_dt > user->dt || (fracpart > 0.001 && fracpart < 0.999))
  {
    PetscPrintf(PETSC_COMM_WORLD,
                "\nThe -ptc_dt value must be less than the -dt value and\n"
                "-ptc_dt must be and integral multiple of -dt\n");
    exit(1);
  }
  ptc_int = intpart;
  if (fracpart >= 0.999)
    ++ptc_int;
  // Calculate most accurate value of ptc_dt.
  ptc_dt = user->dt / ptc_int;

  // Calculate time steps needed to integrate from tis to tie
  ptc_tn = (user->tie - user->tis) * user->dt / ptc_dt;
  PetscPrintf(PETSC_COMM_WORLD, "Time steps to integrate: %i\n", ptc_tn);

  // Change dt from time between each timestep to the time between writting
  // each solution file.
  user->tsdt = user->dt * user->ts;

  // Number of particle iterations between the iterations of the results files.
  ptc_ftn = ptc_int * user->ts;

  // Want to output timing so continue on just don't do important stuff.
  int error = 0;

  // Valid ptc_rstart value is OK.
  if (user->ptc_rstart != -1 && user->ptc_rstart <= user->tis)
  {
    PetscPrintf(PETSC_COMM_WORLD,
            " ptc_rstart value (%d) can't be less than or equal to tis (%d)\n",
              user->ptc_rstart, user->tis);
    error = 1;
  }

  if (!error)
  {
    // Read ib surfaces
    if (user->useIBs)
      new Surfaces(user, user->surfs);
  }
  
  std::vector<Particle> Ptc;

  if (!error)
  {
    PetscPrintf(PETSC_COMM_WORLD, "=============== ParticleInitial\n");
    int ierr = ParticleInitial(user, Ptc);
    if (ierr)
      error = ierr;
  }

  if (!error)
  {
    PetscPrintf(PETSC_COMM_WORLD, "=============== tracking\n");
    int ierr = tracking(user, Ptc);
    if (ierr)
      error = ierr;
  }
  
  OutputGlobalErrors(user);
  
  auto elapTime = user->globTimer.getElapsedCPUTimeStr();
  PetscPrintf(PETSC_COMM_WORLD, "Elapsed time at end of run: %s\n",
              elapTime.c_str());

  PetscFinalize();


  fclose(user->rkdebug);


  return error;
}
