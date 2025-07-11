#include "StringStuff.h"

void splitString(std::string const & s, std::vector < std::string > & split)
{
  size_t currPtr = 0;
  split.resize(0);
  
  while (true)
  {
    size_t startPtr = s.find_first_not_of(" \t", currPtr);
    if (startPtr == std::string::npos)
      break;
    size_t endPtr = s.find_first_of(" \t", startPtr);
    split.push_back(s.substr(startPtr, endPtr-startPtr));
    currPtr = endPtr;
  }
}


std::vector < std::string > splitString(std::string const & s)
{
  std::vector < std::string > split;
  splitString(s, split);
  return split;
}

void splitString(std::string const & s, std::string & split);
std::vector < std::string > splitString(std::string const & s);
