#include "GeomException.h"

namespace Geom
{
  Exception::Exception()
    : CommonException()
  {

  }


  Exception::Exception(std::string const &description)
    : CommonException(description)
  {

  }
}
