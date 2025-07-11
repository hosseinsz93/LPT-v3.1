#ifndef _GEOM_EXCEPTION_H_
#define _GEOM_EXCEPTION_H_


#include <exception>


namespace Geom
{
  class Exception : public std::exception
  {
    public:
      Exception(std::string const & throwClass, bool log) throw()
        : std::exception()
      {
        _mess = throwClass + " thrown.";
      }

      Exception(char const * mess) throw()
        : std::exception()
        , _mess(mess)
      {
      }

      Exception(std::string const & mess) throw()
        : std::exception()
        , _mess(mess)
      {
      }

      ~Exception() throw() { }

      void setDescription(std::string const & desc)
      {
        _mess = desc;
      }

      char const * what() const throw() { return _mess.c_str(); }

    private:
      std::string _mess;
  };
}

#endif
