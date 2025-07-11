#include "GeomId.h"

#include <limits>
#include <iostream>

namespace Geom
{

  Id const Id::min(0);
  Id const Id::firstValid(1);
  Id const Id::max(std::numeric_limits<Id::SignedType>::max());

  std::ostream &operator<< (std::ostream &s, Id const &id)
  {
    if (!id.isValid())
      s << "(invalid)";
    else {
      if (!id.getSign()) s << '-';
      s << id.getValue();
    }
    return s;
  }

  std::istream &operator>> (std::istream &s, Id &id)
  {
    Id::SignedType n;
    s >> n;
    if (!s.fail())
      id = Id(n);
    else
      id = Id::max;
    return s;
  }
}
