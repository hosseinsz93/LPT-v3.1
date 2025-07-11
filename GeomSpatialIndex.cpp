#include <memory>

#include "GeomSpatialIndex.h"
#include "GeomSpatialIndexSource.h"

namespace Geom
{
  SpatialIndex::SpatialIndex(std::shared_ptr<SpatialIndexSource> const source)
    : _source(source->clone())
  {

  }

  SpatialIndex::~SpatialIndex()
  {
  }

}
