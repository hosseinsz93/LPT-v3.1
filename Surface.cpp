#include "Surface.h"

#include <algorithm>
#include <stdio.h>

#include "petsc.h"
#include "GeomRtree.h"
#include "GeomSubset.h"
#include "SpatialIndexBB.h"
#include "SpatialIndexSurf.h"
#include "SpatialIndexSurfs.h"
#include "StringStuff.h"
#include "variables.h"


Surface::Surface(std::string _filename)
  : filename(_filename)
{
  std::fstream fin(filename, std::ios_base::in);
  if (fin.fail())
  {
    PetscPrintf(PETSC_COMM_WORLD, "Can't open file %s\n", filename.c_str());
    exit(1);
  }

  std::string buff;
      
  // The first record contains a list of variables available.  There MUST
  // be at least 3 for the x, y and z coordinates of the mesh.  Count the
  // number of variables now.  The number of ',' +1 is the number of
  // variables.
  getline(fin, buff);
  if (fin.bad() || fin.eof()) fileError();
  remCargRtn(buff);
  int numVar = 1;
  std::size_t pos = 0;
  while (true)
  {
    pos = buff.find_first_of(',', pos);
    if (pos == std::string::npos)
      break;
    ++pos;
    ++numVar;
  }
  if (numVar < 3)
    fileError();

  // The second record contains the number of vertices and cells.
  getline(fin, buff);
  if (fin.bad() || fin.eof()) fileError();
  remCargRtn(buff);
  int cnt = sscanf(buff.c_str(), "ZONE T='TRIANGLES', N=%d, E=%d,",
                   &nverts, &ntris);
  if (cnt != 2) fileError();

  // Allocate space for vectors
  verts.resize(nverts);
  tris.resize(ntris);

  // Read x, y and z coordinates
  readCoor(fin, X);
  readCoor(fin, Y);
  readCoor(fin, Z);

  // Read normals but don't save
  if (numVar >= 4)
    readNorm(fin);
  if (numVar >= 5)
    readNorm(fin);
  if (numVar >= 6)
    readNorm(fin);
  if (numVar >= 7)
    readNorm(fin);
  if (numVar >= 8)
    readNorm(fin);
  if (numVar >= 9)
    readNorm(fin);
  if (numVar >= 10)
    readNorm(fin);
  if (numVar >= 11)
    readNorm(fin);
  if (numVar >= 12)
    readNorm(fin);

  // Read triangles
  readTris(fin);

  for (auto const & v : verts)
    bb.expand(v);

  // outputStar();
}

void Surface::createRtree(std::shared_ptr<Surface> & surf)
{
  // Generate rtree for this ib.
  new SpatialIndexSurf(surf, spatialIndexSurf);
  rtreeSurf.reset(new Geom::Rtree(spatialIndexSurf));

  Geom::Subset set(Geom::Subset::Vector);
  set.resize(Geom::Id(ntris+1));
  for (int i=0; i<ntris; ++i)
    set.add(Geom::Id(i+1));
  rtreeSurf->addSubset(set);
}

void Surface::fileError()
{
  PetscPrintf(PETSC_COMM_WORLD, "Error reading file %s\n",
              filename.c_str());
  exit(1);
}

void Surface::remCargRtn(std::string & buff)
{
  auto result = buff.find_first_of("\r");
  if (result != std::string::npos)
    buff = buff.substr(0, result);
  return;
}
    
void Surface::readCoor(std::fstream & fin, CoorLabel coorIdx)
{
  int cnt = 0;
  std::string buff;
  std::vector < std::string > split;

  while (cnt < nverts)
  {
    getline(fin, buff);
    if (fin.bad() || fin.eof()) fileError();
    remCargRtn(buff);
    splitString(buff, split);
    for (auto const & l : split)
    {
      int c = sscanf(l.c_str(), "%lf", &verts[cnt++](coorIdx));
      if (c != 1)
        fileError();
          
      if (cnt >= nverts)
        break;
    }
  }
}

void Surface::readNorm(std::fstream & fin)
{
  int cnt = 0;
  std::string buff;
  std::vector < std::string > split;

  while (cnt < ntris)
  {
    getline(fin, buff);
    if (fin.bad() || fin.eof()) fileError();
    remCargRtn(buff);
    splitString(buff, split);
    for (auto const & l: split)
      if (cnt++ >= ntris)
        break;
  }
}

void Surface::readTris(std::fstream & fin)
{
  int cnt = 0;
  std::string buff;
  std::vector < std::string > split;

  while (cnt < ntris)
  {
    getline(fin, buff);
    if (fin.bad()) fileError();
    remCargRtn(buff);

    // Empty line
    size_t startPtr = buff.find_first_not_of(" \t");
    if (startPtr == std::string::npos)
      continue;

    splitString(buff, split);
    for (int i=0; i<3; ++i)
    {
      int c = sscanf(split[i].c_str(), "%d", &tris[cnt][i]);
      if (c != 1) fileError();
      --tris[cnt][i];
    }
    ++cnt;
  }
}

void Surface::outputStar() const
{
  char buff[80];

  auto f = fopen("star.vrt", "w");
  for (auto const & v : verts)
  {
    sprintf(buff, "v,%d,%f,%f,%f\n", (int)(&v - &*verts.begin() + 1),
            v(0), v(1), v(2));
    fwrite(buff, 1, strlen(buff), f);
  }
  fclose(f);

  f = fopen("star.cel", "w");
  for (auto const & t : tris)
  {
    sprintf(buff, "c,%d,%d,%d\n", t[0]+1, t[1]+1, t[2]+1);
    fwrite(buff, 1, strlen(buff), f);
  }
  fclose(f);
}

std::vector < Vector3d > Surface::getTriCoor(int tId) const
{
  auto const & t = gettri(tId);
  return std::vector < Vector3d > ({ verts[t[0]], verts[t[1]], verts[t[2]] });
}



Surfaces::Surfaces(std::shared_ptr< UserCtx > & _user,
                   std::shared_ptr< Surfaces > & _self)
  : user(_user)
{
  self.reset(this);
  _self = self;
  
  // If the file surfaces.dat exists, open it and read all the surface file
  // named in the file.  It reads the surface output from vfs.
  std::fstream fin("surfaces.dat", std::ios_base::in);
  if (fin.fail())
  {
    PetscPrintf(PETSC_COMM_WORLD, "No surfaces.dat file found.\n");
    exit(1);
  }

  std::string buff;
  while (true)
  {
    getline(fin, buff);
    if (fin.bad())
    {
      PetscPrintf(PETSC_COMM_WORLD, "Error reading file surfaces.dat.\n");
      exit(1);
    }
    if (fin.eof())
      break;

    Surface::remCargRtn(buff);
    if (buff[0] == '#')
      continue;
    
    std::shared_ptr<Surface> surf(new Surface(buff));
    surf->createRtree(surf);
    surfs.push_back(surf);
  }

  createGlobalRtree();

  createSurfBBRtree();
}

void Surfaces::createGlobalRtree()
{
  vIdxs.resize(surfs.size()+1);
  vIdxs[0] = 0;
  
  tIdxs.resize(surfs.size()+1);
  tIdxs[0] = 0;
  
  for (unsigned s=0; s<surfs.size(); ++s)
  {
    vIdxs[s+1] = vIdxs[s] + surfs[s]->vsize();
    tIdxs[s+1] = tIdxs[s] + surfs[s]->tsize();
  }

  // Generate rtree for triangule faces
  new SpatialIndexSurfs(user, user->spatialIndexSurfs);
  user->rtreeSurfs.reset(new Geom::Rtree(user->spatialIndexSurfs));

  Geom::Subset set(Geom::Subset::Vector);
  set.resize(Geom::Id(tIdxs[tIdxs.size()-1]+1));
  for (int i=0; i<tIdxs[tIdxs.size()-1]; ++i)
    set.add(Geom::Id(i+1));
  user->rtreeSurfs->addSubset(set);
}

void Surfaces::createSurfBBRtree()
{
  new SpatialIndexBB(user->surfs, spatialIndexBB);
  rtreeBB.reset(new Geom::Rtree(spatialIndexBB));

  Geom::Subset set(Geom::Subset::Vector);
  set.resize(Geom::Id(surfs.size()+1));
  for (unsigned i=0; i<surfs.size(); ++i)
    set.add(Geom::Id(i+1));
  rtreeBB->addSubset(set);
}

std::tuple<int,int> Surfaces::convTId(Geom::Id const &id) const
{
  int zId = id.getValue() - 1;
      
  if (tIdxs.size() == 1)
    return std::make_tuple(0, zId);

  
  if (tIdxs.size() < 8)
  {
    unsigned sId;
    for (sId=1; sId<tIdxs.size(); ++sId)
      if (zId < tIdxs[sId])
        break;
    return std::make_tuple(sId-1, zId-tIdxs[sId-1]);
  }

  // 1 added to zId to get the correct surface ib id at the
  // end extreme of each id.
  auto tmp = std::lower_bound(tIdxs.begin(), tIdxs.end(), zId+1);
  int sId = tmp - tIdxs.begin();
  return std::make_tuple(sId-1, zId-tIdxs[sId-1]);
}


std::vector < Vector3d > Surfaces::getTriCoor(Geom::Id const &id) const
{
  int sId=0, zId=0;
  std::tie(sId, zId) = convTId(id);
  return getSurf(sId)->getTriCoor(Geom::Id(zId+1));
}
