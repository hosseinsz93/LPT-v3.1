#include "GeomCompressionData.h"
#include "GeomVectorStorage.h"

namespace Geom
{
  CompressionData::CompressionData()
    : _maxUncompressedId(Id::min),
      _maxCompressedId(Id::min)
  {}

  CompressionData::~CompressionData()
  {}

  void CompressionData::clear()
  {
    std::vector<Id>().swap(_uncompressedToCompressed);

    _maxUncompressedId = Id::min;
    _maxCompressedId = Id::min;
  }

  void CompressionData::initialize(Id maxUncompressedId,
				   Id maxCompressedId)
  {
    _maxUncompressedId = maxUncompressedId;
    _maxCompressedId = maxCompressedId;
    
    _uncompressedToCompressed.resize(maxUncompressedId.toUnsigned() + 1);

    std::fill(_uncompressedToCompressed.begin(),
	       _uncompressedToCompressed.end(),
	       Id::max);
  }

  void CompressionData::buildInverseMap(VectorStorage<Id> &inverse) const
  {
    inverse.resize(_maxCompressedId);
    for (unsigned i = 0; i<_uncompressedToCompressed.size(); ++i)
      if (_uncompressedToCompressed[i] != Id::max)
	// do range checking here incase map has a bogus value...
	inverse.at(_uncompressedToCompressed[i]) = Id(i);
  }
}
