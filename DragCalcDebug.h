#ifndef _DRAG_CALC_DEBUG_H_
#define _DRAG_CALC_DEBUG_H_

#include "variables.h"

#include "GeomVector.h"

#include <stdio.h>


class DragCalcDebug
{
  public:
    DragCalcDebug(int _fileIdx=-1);
    ~DragCalcDebug();
    void swapFile(int _fileIdx=-1);
    void close();
    void add(std::shared_ptr< UserCtx > & user,
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
             int partStat);

  private:
    int fileIdx;
    FILE * fd;
};

#endif
