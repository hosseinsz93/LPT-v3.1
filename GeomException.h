#ifndef _GEOM_EXCEPTION_H
#define _GEOM_EXCEPTION_H

#include "GeomCommonException.h"

namespace Geom
{
  class Id;

  class Exception : public CommonException
  {
  public:
    Exception();
    Exception(std::string const &description);
    ~Exception() throw() {}
  };
}

#endif
