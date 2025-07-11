#include "variables.h"

#include <stdio.h>

#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "Id.h"
#include "SpatialIndexGrid.h"
#include "GeomRtree.h"
#include "GeomSubset.h"
#include "GeomXTol.h"
#include "StringStuff.h"

extern PetscInt ptc_tn;
extern PetscInt ptc_ftn;
extern PetscReal cl;

PetscErrorCode ParticleInitialRestart(std::shared_ptr< UserCtx > & user,
                                      std::vector<Particle> & Ptc);


void genRtree(std::shared_ptr< UserCtx > & user)
{
  PetscPrintf(PETSC_COMM_WORLD, "DEBUG starting grid genRtree\n");
  Util::Timer timer(Util::Timer::Wall);

  // Generate rtree for grid;
  new SpatialIndexGrid(user, user->spatialIndexGrid);
  user->rtreeGrid.reset(new Geom::Rtree(user->spatialIndexGrid));

  Geom::Subset set(Geom::Subset::Vector);
  set.resize(Id(user->IM, user->JM, user->KM));
  PetscInt xs, xe, ys, ye, zs, ze;
  user->spatialIndexGrid->indexRanges(xs, xe, ys, ye, zs, ze);
  for (int k=zs; k<ze-2; ++k)
  for (int j=ys; j<ye-2; ++j)
  for (int i=xs; i<xe-2; ++i)
  {
    set.add((Geom::Id) Id(i, j, k));
  }
  user->rtreeGrid->addSubset(set);

  std::stringstream mess;
  mess << timer;
  PetscPrintf(PETSC_COMM_WORLD, "%s\n", mess.str().c_str());
  PetscPrintf(PETSC_COMM_WORLD, "DEBUG ending grid genRtree\n");
}


PetscErrorCode Initial(FILE *fd, std::shared_ptr< UserCtx > & user,
                       Util::Timer & timer)
{
  PetscPrintf(PETSC_COMM_WORLD, "block_number %d\n", block_number);
  
  PetscInt bi;
  PetscReal *gc;
  PetscInt i, j, k;
  //  PetscReal cl = 19.5;

  // bi MUST always be zero because shared_ptr's are used.  You would
  // need a vector of shared_ptr's for this to work.
  //for (bi=0; bi<block_number; bi++)
  bi = 0;
  {
    if (binary_input)
    {
      fread(&user->IM, sizeof(int), 1, fd);
      fread(&user->JM, sizeof(int), 1, fd);
      fread(&user->KM, sizeof(int), 1, fd);
    }
    else
      fscanf(fd, "%i %i %i\n", &user->IM, &user->JM, &user->KM);
    
    PetscInt IM, JM, KM;
    IM = user->IM; JM = user->JM; KM = user->KM;
    PetscPrintf(PETSC_COMM_WORLD, "block: %d  i,j,k: %d,%d,%d\n",
                bi, IM, JM, KM);

    DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_BOX,
	       user->IM+1, user->JM+1, user->KM+1,
	       PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, 2, PETSC_NULL,
	       PETSC_NULL, PETSC_NULL, &user->da);
    DASetUniformCoordinates(user->da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    DAGetCoordinateDA(user->da, &user->fda);
    DAGetLocalInfo(user->da, &user->info);

    PetscMalloc(3*(IM*JM*KM)*sizeof(PetscReal), &gc);

    for (k=0; k<KM; k++) {
    for (j=0; j<JM; j++) {
	for (i=0; i<IM; i++) {
          if (binary_input)
            fread(gc + (k*(JM*IM) + j * IM + i)*3, sizeof(double), 1, fd);
          else
            fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3);
	}
    }
    }
      
    for (k=0; k<KM; k++) {
      for (j=0; j<JM; j++) {
	for (i=0; i<IM; i++) {
          if (binary_input)
            fread(gc + (k*(JM*IM) + j * IM + i)*3 + 1, sizeof(double), 1, fd);
          else
            fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3 + 1);
	}
      }
    }

    for (k=0; k<KM; k++) {
      for (j=0; j<JM; j++) {
	for (i=0; i<IM; i++) {
          if (binary_input)
            fread(gc + (k*(JM*IM) + j * IM + i)*3 + 2, sizeof(double), 1, fd);
          else
            fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3 + 2);
	}
      }
    }

    DALocalInfo	info = user->info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;

    Vec Coor;
    Cmpnts ***coor;

    DAGetGhostedCoordinates(user->da, &Coor);
    DAVecGetArray(user->fda, Coor, &coor);
    for (k=zs; k<ze; k++)
    for (j=ys; j<ye; j++)
    for (i=xs; i<xe; i++)
      if (k<KM && j<JM && i<IM) 
      {
        coor[k][j][i].x = *(gc + (k * (IM*JM) + j * IM + i) * 3  )/cl;
        coor[k][j][i].y = *(gc + (k * (IM*JM) + j * IM + i) * 3+1)/cl;
        coor[k][j][i].z = *(gc + (k * (IM*JM) + j * IM + i) * 3+2)/cl;
      }

    PetscFree(gc);

    DACreateGlobalVector(user->fda, &user->Ucat);
    VecDuplicate(user->Ucat, &user->Ucat_c);
    VecDuplicate(user->Ucat, &user->Cent);
    VecDuplicate(user->Ucat, &user->Ucat_o);
    
    DACreateGlobalVector(user->da, &user->Rad);

    Cmpnts ***cent;

    DAVecGetArray(user->fda, user->Cent, &cent);

    // calculate cell centroid
    for (k=zs; k<mz-2; k++)
    for (j=ys; j<my-2; j++)
    for (i=xs; i<mx-2; i++)
      cent[k][j][i] = CENT(coor, k, j, i);
   
    DAVecRestoreArray(user->fda, user->Cent, &cent);
    DAVecRestoreArray(user->fda, Coor, &coor);
  }
  fclose(fd);

  std::stringstream mess;
  mess << timer;
  PetscPrintf(PETSC_COMM_WORLD, "%s\n", mess.str().c_str());
  PetscPrintf(PETSC_COMM_WORLD, "DEBUG finished reading grid.dat\n");

  return 0;
}


PetscErrorCode ParticleInitial(std::shared_ptr< UserCtx > & user,
                               std::vector<Particle> & Ptc)
{
  if (user->initPart == 1)
  {
  /* For the time being, particles are released at the center of each grid cells at a constant k plane */
  /* PetscInt num = -1, k_slides=20; */

  DALocalInfo   info = user->info;
  PetscInt      xs = info.xs, xe = info.xs + info.xm;
  PetscInt      ys = info.ys, ye = info.ys + info.ym;
  PetscInt      zs = info.zs, ze = info.zs + info.zm; //Ali added 4/15/15

//Ali -> Maximum number of Particles must be equal to (xe-2)*(ye-2): e.g. for a BackGround grid of 101 X 53 the maximum particle number, ptc_num, should be = 100 X 52
//  PetscPrintf(PETSC_COMM_WORLD, "k_slides, ptc_num, total_point %d %d %d \n",k_slides,user->ptc_num, pointstotal);
  PetscPrintf(PETSC_COMM_WORLD, "xs xe xm %d %d %d \n",xs,xe,info.xm);
  PetscPrintf(PETSC_COMM_WORLD, "ys ye ym %d %d %d \n",ys,ye,info.ym);
  PetscPrintf(PETSC_COMM_WORLD, "zs ze zm %d %d %d \n",zs,ze,info.zm); //Ali added 4/15/15

  double x, y, z;
  double u, v, w;
  double den, dia, stime;

  std::vector <std::string> tmp_particle;

  std::ifstream pin("ParticleInitial.dat", std::ifstream::in);
  if (!pin.is_open())
  {
    PetscPrintf(PETSC_COMM_WORLD, "Can't open file ParticleInitial.dat.\n");
    return 1;
  }

  std::string line;
  while (std::getline(pin, line))
  {
    if (line[0] == '#')
      continue;
    int cnt = sscanf(line.c_str(), "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
           &x, &y, &z, &u, &v, &w, &stime, &dia, &den);
    if (!cnt)
      continue;
    tmp_particle.push_back(line);
  }

  pin.close();


  user->ptc_num = tmp_particle.size();
  Ptc.resize(user->ptc_num);

  auto fd = fopen("Particle_Init.dat", "w");
  PetscFPrintf(PETSC_COMM_WORLD, fd,
               "Variables = X,\tY,\tZ,\tU,\tV,\tW,\tSTIME,\tETIME,\tDIA,\tDEN,\tOUT\n");
  PetscFPrintf(PETSC_COMM_WORLD, fd,
               "ZONE I = %i,\tJ = %i, DATAPACKING = POINT\n", 1, user->ptc_num);

  Cmpnts ***cent;
  DAVecGetArray(user->fda, user->Cent, &cent);

  for (auto const & l : tmp_particle)
  {
    auto idx = &l - &tmp_particle.front();

    auto split = splitString(l);

    int cnt = sscanf(l.c_str(), "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
           &x, &y, &z, &u, &v, &w, &stime, &dia, &den);


#if 0
    switch (cnt)
    {
      case  1:
        PetscPrintf(PETSC_COMM_WORLD, "DEBUG ptc values read: x: %f\n", x);
        break;

      case  2:
        PetscPrintf(PETSC_COMM_WORLD,
                    "DEBUG ptc values read: x,y: %f,%f\n", x, y);
        break;

      case  3:
        PetscPrintf(PETSC_COMM_WORLD,
                    "DEBUG ptc values read: x,y,z: %f,%f,%f\n", x, y, z);
        break;

      case  4:
        PetscPrintf(PETSC_COMM_WORLD,
                    "DEBUG ptc values read: x,y,z: %f,%f,%f  "
                    "u: %f\n", x, y, z, u);
        break;

      case  5:
        PetscPrintf(PETSC_COMM_WORLD,
                    "DEBUG ptc values read: x,y,z: %f,%f,%f  "
                    "u,v: %f,%f\n", x, y, z, u, v);
        break;

      case  6:
        PetscPrintf(PETSC_COMM_WORLD,
                    "DEBUG ptc values read: x,y,z: %f,%f,%f  "
                    "u,v,w: %f,%f,%f\n",
                    x, y, z, u, v, w);
        break;

      case  7:
        PetscPrintf(PETSC_COMM_WORLD,
                    "DEBUG ptc values read: x,y,z: %f,%f,%f  "
                    "u,v,w: %f,%f,%f  stime: %f\n",
                    x, y, z, u, v, w, stime);
        break;

      case  8:
        PetscPrintf(PETSC_COMM_WORLD,
                    "DEBUG ptc values read: x,y,z: %f,%f,%f  "
                    "u,v,w: %f,%f,%f  stime,dia: %f,%e\n",
                    x, y, z, u, v, w, stime, dia);
        break;

      case  9:
        PetscPrintf(PETSC_COMM_WORLD,
                    "DEBUG ptc values read: x,y,z: %f,%f,%f  "
                    "u,v,w: %f,%f,%f  stime,dia,den: %f,%e,%f\n",
                    x, y, z, u, v, w, stime, dia, den);
        break;
    }
#endif
    

    // Check that the first 3 characters are integers in which case
    // it is the i,j,k values that were input.
    if (cnt >= 3)
    {
      int ijk[3] = { 0, 0, 0 };
      bool areInt = true;

      for (int i=0; i<3; ++i)
      {
        char * endptr;
        char const * str = split[i].c_str();
        int value = strtol(str, &endptr, 10);
        if (*endptr != 0)
        {
          areInt = false;
          break;
        }
      
        switch (i)
        {
          case 0:
            ijk[i] = Ptc[idx].loc.i = value;
            break;

          case 1:
            ijk[i] = Ptc[idx].loc.j = value;
            break;

          case 2:
            ijk[i] = Ptc[idx].loc.k = value;
            break;
        }
      }

      if (areInt)
      {
        auto const & c = cent[ijk[2]] [ijk[1]] [ijk[0]];
        x = c.x;
        y = c.y;
        z = c.z;
      }
    }
      
    switch (cnt)
    {
      case 9:
        Ptc[idx].DDenPt = den;
      case 8:
        Ptc[idx].DDiaPt = dia;
      case 7:
        Ptc[idx].stime  = stime;
      case 6:
        Ptc[idx].Vp.set(2, w);
      case 5:
        Ptc[idx].Vp.set(1, v);
      case 4:
        Ptc[idx].Vp.set(0, u);
      case 3:
        Ptc[idx].coor.z = z;
      case 2:
        Ptc[idx].coor.y = y;
      case 1:
        Ptc[idx].coor.x = x;
    }

    if (cnt < 8)
      Ptc[idx].DDiaPt = user->DDiaPt;
    if (cnt < 9)
      Ptc[idx].DDenPt = user->DRoPt;

    Ptc[idx].etime = Ptc[idx].stime;

    PetscFPrintf(PETSC_COMM_WORLD, fd,
                 "%d\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%d\n",
                 idx, Ptc[idx].coor.x, Ptc[idx].coor.y, Ptc[idx].coor.z,
                 Ptc[idx].Vp(0), Ptc[idx].Vp(1), Ptc[idx].Vp(2), Ptc[idx].stime,
                 Ptc[idx].etime, Ptc[idx].DDiaPt, Ptc[idx].DDenPt, Ptc[idx].partStat);
  }

  DAVecRestoreArray(user->fda, user->Cent, &cent);

  fclose(fd);
  }

  if (user->ptc_rstart > 0)
  {
    auto ierr = ParticleInitialRestart(user, Ptc);
    if (ierr)
      return ierr;
  }

  PetscInt tmptis = (user->ptc_rstart > 0) ? user->ptc_rstart : user->tis;
  return VelocityReadin(user, tmptis, 0);
}


PetscErrorCode ParticleInitialRestart(std::shared_ptr< UserCtx > & user,
                                      std::vector<Particle> & Ptc)
{
  char fbuff[128] = "";
  sprintf(fbuff, "Part.Traj.%06d.dat", (int)user->ptc_rstart);

  std::ifstream pin(fbuff, std::ifstream::in);
  if (!pin.is_open())
  {
    PetscPrintf(PETSC_COMM_WORLD, "Can't open restart file %s\n", fbuff);
    return 1;
  }

  bool timeFound = false;
  double steptime = user->ptc_rstart * user->dt;
  std::string line;
  while (std::getline(pin, line))
  {
    if (line[0] == '#')
      continue;

    if (line.substr(0,4) == "time")
      continue;
        
    int partStat, pId;
    double x, y, z, u, v, w, stime, etime, dia, den, time;
    int cnt = sscanf(line.c_str(),
                     "%lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",
                     &time, &pId, &x, &y, &z, &u, &v, &w, &stime, &etime, &dia,
                     &den, &partStat);
    if (cnt != 13)
    {
      PetscPrintf(PETSC_COMM_WORLD, "Error reading file %s\n", fbuff);
      return 1;
    }

    // There are many time entries stored in this file.  Only use the one
    // for the time step.
    if (time != GeomX::Tol(steptime))
      continue;
    
    Ptc[pId].coor.x = x;
    Ptc[pId].coor.y = y;
    Ptc[pId].coor.z = z;
    Ptc[pId].Vp.set(u, v, w);
    
    // std::vector < Cmpnts > traj;
    // PetscReal dt; /* time step used to integrate the particle trajectory */
    // PetscInt next_time;
    Ptc[pId].firstVel = true;

    Ptc[pId].stime = stime;
    Ptc[pId].etime = etime;
    Ptc[pId].DDiaPt = dia;
    Ptc[pId].DDenPt = den;

    Ptc[pId].dist = 0.0;
    
    Ptc[pId].partStat = static_cast<ParticleStatus>(partStat);

    // These values must be set correctly before going into ParticleLocation
    // subroutine.  The show the status of the geometry lookup to which finds
    // the i,j,k value for the particles current location.
    // Allowable values are defined in enum GeomLookupStatus.
    // GEOM_NOT_DONE,
    // GEOM_SUCCESSFUL, The closes cell has been found and the i,j,k is set.
    // GEOM_FAILED      Geom Lookup failed to find a cell close enough.
    Ptc[pId].geomLookupStatus = GEOM_NOT_DONE;

    timeFound = true;
  }
  
  pin.close();

  if (! timeFound)
  {
    PetscPrintf(PETSC_COMM_WORLD,
                "No entries found for the restart time in file %s.\n"
                "Try using a previous time step.\n", fbuff);
    return 1;
  }
  
  return 0;
}


