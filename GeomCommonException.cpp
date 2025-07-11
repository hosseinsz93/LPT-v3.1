#include "GeomCommonException.h"

#include <sstream>

namespace Geom
{
  CommonException::CommonException()
    : _description("(undefined exception)")
  {

  }

  CommonException::CommonException(std::string const &description)
    : _description(description)
  {

  }


  char const *CommonException::what() const throw()
  {
    return _description.c_str(); 
  }

  void CommonException::setDescription(std::string const &description)
  {
    _description = description;
  }


  InvalidIndex::InvalidIndex(Id const &n)
  {
    std::ostringstream str;
    str << "Invalid Index ";
    str << n;
    setDescription(str.str());
  }

  
  MasterSerialException::MasterSerialException(std::string const &description)
    : CommonException(description)
  {}
}
