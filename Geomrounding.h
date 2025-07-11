#ifndef _MKGEOMROUNDING_H
#define _MKGEOMROUNDING_H

// fpu control

#if (defined(_AIX)      || defined(_IBMR2)    || defined(__IBMR2__) || \
     defined(__alpha)   || defined(__alpha__) ||                       \
     defined(hpux)      || defined(__hpux)    || defined(__hpux__)  || \
     defined(__ia64)    || defined(__ia64__)  ||                       \
     defined(sgi)       || defined(__sgi)     || defined(__sgi__)   || \
     defined(sun)       || defined(__sun)     || defined(__sun__)   || \
    (defined(__PPC64__) && defined(__linux__)))
#define FPU_ROUND_DOUBLE
#define FPU_RESTORE

#elif (defined(WIN32) || defined(WIN64))
// On Windows x64 architecture, changing the floating point precision is not supported.
//# if defined(WIN64)
#define FPU_ROUND_DOUBLE
#define FPU_RESTORE
//# else
//#include <float.h>
//static unsigned int fpu_init;
//#define FPU_ROUND_DOUBLE {fpu_init = _controlfp (0, 0); _controlfp (_PC_53, _MCW_PC);}
//#define FPU_RESTORE {_controlfp (fpu_init, 0xfffff);}
//# endif

#else
/* FPU control. We MUST have only double precision (not extended precision) */
#include <fpu_control.h>
static fpu_control_t fpu_round_double =
      (_FPU_DEFAULT & ~ _FPU_EXTENDED)|_FPU_DOUBLE;
static fpu_control_t fpu_init;
#define FPU_ROUND_DOUBLE  {_FPU_GETCW(fpu_init);_FPU_SETCW(fpu_round_double);}
#define FPU_RESTORE       {_FPU_SETCW(fpu_init);}
#define USE_PREDICATES_INIT
#endif

// exact control

#if (defined(hpux) || defined(__hpux) || defined(__hpux__))
#define INEXACT volatile
#else
#define INEXACT
#endif

#endif
