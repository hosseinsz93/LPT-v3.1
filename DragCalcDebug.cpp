#include "DragCalcDebug.h"


DragCalcDebug::DragCalcDebug(int _fileIdx)
  : fileIdx(0)
  , fd(NULL)
{
  if (_fileIdx != -1)
    swapFile(_fileIdx);
}


DragCalcDebug::~DragCalcDebug()
{
  close();
}


void DragCalcDebug::swapFile(int _fileIdx)
{
  close();

  char fname[128];
  if (_fileIdx != -1)
    fileIdx = _fileIdx;
  else
    ++fileIdx;
  sprintf(fname, "dragCalcDebug%06d.dat", fileIdx);
  fd = fopen(fname, "w");
  if (!fd)
  {
    PetscPrintf(PETSC_COMM_WORLD, "Couldn't open file %s\n", fname);
    exit(1);
  }

  fprintf(fd,
    /* int const pId */          "pId"
    /* double const time */      "\ttime"
    /* Vector3d const P0 */      "\tP0x\tP0y\tP0z"
    /* Location const loc0 */    "\tloc0.i\tloc0.j\tloc0.k"
    /* Vector3d const Vf0 */     "\tVf0x\tVf0y\tVf0z\tVf0.mag"
    /* Vector3d const Vp0 */     "\tVp0x\tVp0y\tVp0z\tVp0.mag"
    /* Vector3d const Vpf0 */    "\tVpf0x\tVpf0y\tVpf0z\tVpf0.mag"
    /* double const DDiaPt */    "\tDDiaPt"
    /* double const Fden */      "\tDDenFd"
    /* double const DDenPt */    "\tDDenPt"
    /* double const Masp */      "\tDMasPt"
    /* double const Fgrav */     "\tDAccGrav"
    /* double const rey */       "\trey"
    /* double const resptime */  "\trespTime"
    /* double const Cd */        "\tCd"
    /* Vector3d const Fdrag1 */  "\tFdrag1.mag"
    /* Vector3d const Fdrag2 */  "\tFdrag2.mag"
    /* Vector3d const Fdrag3 */  "\tFdrag3.mag"
    /* Vector3d const Fdrag */   "\tFdragx\tFdragy\tFdragz\tFdrag.mag"
    /* Vector3d const & DFtot */ "\tDFtotx\tDFtoty\tDFtotz\tDFtot.mag"
    /* Vector3d const & Dacc */  "\tDaccx\tDaccy\tDaccz\tDacc.mag"
    /* Vector3d const & DVel */  "\tDVelx\tDVely\tDVelz\tDVel.mag"
    /* double const Vt */        "\tVt"
    /* compare Vpf1 to Vt */     "\tVt/Vpf1%"
    /* Vector3d const Vf1 */     "\tVf1x\tVf1y\tVf1z\tVf1.mag"
    /* double const Vp1 */       "\tVp1x\tVp1y\tVp1z\tVp1.mag"
    /* double const Vpf1 */      "\tVpf1x\tVpf1y\tVpf1z\tVpf1.mag"
    /* Vector3d const P1 */      "\tP1x\tP1y\tP1z"
    /* Location const loc1 */    "\tloc1.i\tloc1.j\tloc1.k"
    /* int const partStat */     "\tpStat"
                                 "\n");
}


void DragCalcDebug::close()
{
  if (fd)
  {
    fclose(fd);
    fd = NULL;
  }
}


void DragCalcDebug::add(std::shared_ptr< UserCtx > & user,
                        int const pId,
                        double const time,
                        Vector3d const & Vf0,
                        Vector3d const & Vp0,
                        Vector3d const & Vpf0,
                        double const DDiaPt,
                        double const DDenFd,
                        double const DDenPt,
                        double const DMasPt,
                        double const DAccGrav,
                        double const rey,
                        double const respTime,
                        double const Cd,
                        Vector3d const & DFdrag1,
                        Vector3d const & DFdrag2,
                        Vector3d const & DFdrag3,
                        Vector3d const & DFdrag,
                        Vector3d const & DFtot,
                        Vector3d const & Dacc,
                        Vector3d const & DVel,
                        double const Vt,
                        Vector3d const & Vf1,
                        Vector3d const & Vp1,
                        Vector3d const & Vpf1,
                        Vector3d const & P0,
                        Location const & loc0,
                        Vector3d const & P1,
                        Location const & loc1,
                        int partStat)



{
  if (user->useDebugDrag < 0 ||
      (user->useDebugDrag > 0 && user->useDebugDrag == pId))
  {
    fprintf(fd,
      /* int const pId */          "%d"
      /* double const time */      "\t%15e"
      /* Vector3d const P0 */      "\t%15e\t%15e\t%15e"
      /* Location const loc0 */    "\t%d\t%d\t%d"
      /* Vector3d const Vf0 */     "\t%15e\t%15e\t%15e\t%15e"
      /* Vector3d const Vp0 */     "\t%15e\t%15e\t%15e\t%15e"
      /* Vector3d const Vpf0 */    "\t%15e\t%15e\t%15e\t%15e"
      /* double const DDiaPt */    "\t%15e"
      /* double const DDenFd */    "\t%15e"
      /* double const DDenPt */    "\t%15e"
      /* double const DMasPt */    "\t%15e"
      /* double const DAccGrav */  "\t%15e"
      /* double const rey */       "\t%15e"
      /* double const Cd */        "\t%15e"
      /* Vector3d const DFdrag1 */ "\t%15e"
      /* Vector3d const DFdrag2 */ "\t%15e"
      /* Vector3d const DFdrag3 */ "\t%15e"
      /* Vector3d const DFdrag */  "\t%15e\t%15e\t%15e\t%15e"
      /* Vector3d const & DFtot */ "\t%15e\t%15e\t%15e\t%15e"
      /* Vector3d const & Dacc */  "\t%15e\t%15e\t%15e\t%15e"
      /* Vector3d const & DVel */  "\t%15e\t%15e\t%15e\t%15e"
      /* double const Vt */        "\t%15e"
      /* compare Vpf1 to Vt */     "\t%f.4"
      /* Vector3d const Vf1 */     "\t%15e\t%15e\t%15e\t%15e"
      /* double const Vp1 */       "\t%15e\t%15e\t%15e\t%15e"
      /* double const Vpf1 */      "\t%15e\t%15e\t%15e\t%15e"
      /* Vector3d const P1 */      "\t%15e\t%15e\t%15e"
      /* Location const loc1 */    "\t%d\t%d\t%d"
      /* int const partStat */     "\t%d"
      "\n",

      /* int const pId */          pId,
      /* double const time */      time,
      /* Vector3d const P0 */      P0(0), P0(1), P0(2),
      /* Location const loc0 */    loc0.i, loc0.j, loc0.k,
      /* Vector3d const Vf0 */     Vf0(0), Vf0(1), Vf0(2), Vf0.mag(),
      /* Vector3d const Vp0 */     Vp0(0), Vp0(1), Vp0(2), Vp0.mag(),
      /* Vector3d const Vpf0 */    Vpf0(0), Vpf0(1), Vpf0(2), Vpf0.mag(),
      /* double const DDiaPt */    DDiaPt,
      /* double const DDenFd */    DDenFd,
      /* double const DDenPt */    DDenPt,
      /* double const DMasPt */    DMasPt,
      /* double const DAccGrav */  DAccGrav,
      /* double const rey */       rey,
      /* double const respTime */  respTime,
      /* double const Cd */        Cd,
      /* Vector3d const DFdrag1 */ DFdrag1.mag(),
      /* Vector3d const DFdrag2 */ DFdrag2.mag(),
      /* Vector3d const DFdrag3 */ DFdrag3.mag(),
      /* Vector3d const DFdrag */  DFdrag(0), DFdrag(1), DFdrag(2),
                                                                  DFdrag.mag(),
      /* Vector3d const & DFtot */ DFtot(0), DFtot(1), DFtot(2), DFtot.mag(),
      /* Vector3d const & Dacc */  Dacc(0), Dacc(1), Dacc(2), Dacc.mag(),
      /* Vector3d const & DVel */  DVel(0), DVel(1), DVel(2), DVel.mag(),
      /* double const Vt */        Vt,
      /* compare Vpf1 to Vt */     ((DVel.mag() / Vt) - 1) * 100,
      /* Vector3d const Vf1 */     Vf1(0), Vf1(1), Vf1(2), Vf1.mag(),
      /* double const Vp1 */       Vp1(0), Vp1(1), Vp1(2), Vp1.mag(),
      /* double const Vpf1 */      Vpf1(0), Vpf1(1), Vpf1(2), Vpf1.mag(),
      /* Vector3d const P1 */      P1(0), P1(1), P1(2),
      /* Location const loc1 */    loc1.i, loc1.j, loc1.k,
      /* int const partStat */     partStat);

    fflush(fd);
  }
}
