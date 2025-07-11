#ifndef _GEOM_COMMON_EXCEPTION
#define _GEOM_COMMON_EXCEPTION

#include <exception>
#include <string>

#include "GeomId.h"

namespace Geom
{
  class CommonException : public std::exception
  {
  public:
    CommonException();
    CommonException(std::string const &description);
    virtual ~CommonException() throw() {}
    virtual char const *what() const throw();
  protected:
    void setDescription(std::string const &description);
  private:
    std::string _description;
  };

  class InvalidIndex : public CommonException
  {
  public:
    InvalidIndex(Id const &id);
    ~InvalidIndex() throw() {}
  };
  
  /**
   * This class is to be used only in code that is executed only by the master
   * or serial processes.  Any thrown objects should be caught before leaving
   * the master-serial block.  Thus, if an object is thrown from a method,
   * the call graph for that method should be traced upward to guarantee that
   * this exception is caught before the master-serial block is exited.
   * Commenting that the method throws this exception in its API would also
   * be helpful in informing consumers.
   */
  class MasterSerialException : public CommonException
  {
  public:
    MasterSerialException(std::string const &description);
    ~MasterSerialException() throw() {}
  private:
    MasterSerialException();
  };
}

#endif
