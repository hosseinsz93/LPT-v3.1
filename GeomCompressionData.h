#ifndef _GEOM_COMPRESSIONDATA
#define _GEOM_COMPRESSIONDATA

#include <vector>

#include "GeomId.h"

namespace Geom
{
  template<class T> class VectorStorage;

  /*! The CompressionData class is used to generate and maintain a mapping
    from old to new ids during the act of compression.  It is populated by
    a call to Subset::compressMaster() and can then be passed to a variety of
    other containers (subsets, id controllers, storage, etc.) so they can
    update their Ids based on the new data.       
  */
  class CompressionData
  {
  public:
    CompressionData();
    ~CompressionData();

    /*! Get the compressed (new) Id corresponding to a given uncompressed (old)
      Id.  If the Id passed in is invalid or does not correspond to an Id in
      the map, an invalid Id will be returned.  
    */
    Id getCompressedId(Id uncompressedId) const {
      if (!uncompressedId.isValid() ||
	  uncompressedId.toUnsigned() >= _uncompressedToCompressed.size())
	return Id::max;
      return _uncompressedToCompressed[uncompressedId.toUnsigned()];
    }

    /*! Get the maximum compressed (new) Id resulting from the compression. */
    Id getMaxCompressedId() const {
      return _maxCompressedId;
    }

    /*! Get the maximum uncompressed (old) Id from before the compression. */
    Id getMaxUncompressedId() const {
      return _maxUncompressedId;
    }

    /*! Clear the CompressionData object. */
    void clear();

    /*! Initialize the CompressionData object. */
    void initialize(Id maxUncompressedId,
		    Id numElements);

    /*! Set the compressed Id corresponding to an uncompressed Id. */
    void setCompressedId(Id uncompressedId, Id compressedId) {
      _uncompressedToCompressed.at(uncompressedId.toUnsigned()) = 
	compressedId;
    }

    /*! Build the inverse map (compressed->uncompressed ids) and store in the
      input VectorStorage. */
    void buildInverseMap(VectorStorage<Id> &inverse) const;

    
  private:
    Id _maxUncompressedId;
    Id _maxCompressedId;

    std::vector<Id> _uncompressedToCompressed;
  };
}

#endif
