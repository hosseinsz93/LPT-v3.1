#ifndef _GEOM_VECTORSTORAGE
#define _GEOM_VECTORSTORAGE

#include <vector>

#include "GeomUtility.h"

#include "GeomId.h"
#include "GeomStorage.h"
#include "GeomException.h"
#include "GeomCompressionData.h"

namespace Geom
{
  /*!
    A simple concrete class of Storage which wraps around an std::vector and
    provides Id-keyed access to elements.  Since it wraps around a std::vector,
    it will resize itself using whatever mechanism std::vector does, i.e. size
    doubling, and will hold all of it's data in a contiguous array.  If size
    doubling is not appropriate and a contiguous array is not required, use
    BlockedVectorStorage instead.
   */
  template < class T > class VectorStorage : public Storage
  {
  public:
    typedef typename std::vector<T> tbltype;
    typedef typename tbltype::value_type value_type;
    typedef typename tbltype::reference reference;
    typedef typename tbltype::const_reference const_reference;
        
    VectorStorage() {}              
    VectorStorage(Id const &id) {
      _data.resize(id.toUnsigned() + 1);
    }
    ~VectorStorage() {}

    VectorStorage(VectorStorage const &rhs)
      : Storage(rhs)
      , _data(rhs._data)
    {}

    VectorStorage &operator=(VectorStorage const &rhs) {
      if (&rhs == this)
	return *this;

      _data = rhs._data;
      return *this;
    }
    
    /*!
      Checks that the index is valid for the size of the array.
      
      This routine checks that the index is valid in term of the size of the
      array.  In other words, it checks that the index is within the bounds
      of the array.
    */
    bool validId (Id const &idx) const {        
      return (idx.isValid() &&
              idx.toUnsigned() <= maxId().toUnsigned());
    }
        

    /*!
      Resize the container to hold id.
    */
    void resizeNoShrink(Id const &id) {
      if (id.toUnsigned() >= _data.size())
	_data.resize(id.toUnsigned() + 1);
    }

    /*!
      Same as above resize, but specify a value for all new elements.
    */
    void resizeNoShrink(Id const &id, T defaultval) {
      if (id.toUnsigned() >= _data.size())
	_data.resize(id.toUnsigned() + 1, defaultval);
    }

    /*!
      Resize the container to hold id.
    */
    void resize(Id const &id) {      
      _data.reserve(id.toUnsigned() + 1);
      _data.resize(id.toUnsigned() + 1);
    }
    
    /*!
      Same as above resize, but specify a value for all new elements.
    */
    void resize(Id const &id, T defaultval) {
      _data.resize(id.toUnsigned() + 1, defaultval);
    }
    

    /*!
      Is the container empty?
    */
    bool empty() const {
      return _data.empty();
    }

    /*!
      Get size of array.  Note the difference between size and maxId.
      If you want to resize another VectorStorage object, use maxId() instead.
     */
    unsigned size() const {
      return _data.size();
    }

    /*!
      Return the maximum valid id in the container based on its
      current size.
    */
    Id maxId() const { 
      if (_data.empty())
        return Id::min;

      return Id(_data.size() - 1); 
    }


    /*!
      Put the information associated with idx into the container,
      first making sure the array is sized to allow this.
    */
    void resizeAndSet (Id const &id, T const &arg) {
      if (id.toUnsigned() >= _data.size())
	_data.resize(id.toUnsigned() + 1);
      _data[id.toUnsigned()] = arg;
      return;
    }

    /*!
      Put the information associated with id into the container,
      without resizing the array.
    */
    void set(Id const &id, T const &arg) {
      _data[id.toUnsigned()] = arg;
    }


    /// Set all values in the array to val.
    void setAll(T const &arg) {
      std::fill(_data.begin(), _data.end(), arg);
    }


    ///  non-const, non-range checked access into the array
    reference operator[] (unsigned const &id) {
      return _data[id];
    }
    reference operator[] (Id const &id) {
      return (*this)[id.toUnsigned()];
    }

    /// const, non-range checked access into the array
    const_reference operator[] (unsigned const &id) const {
      return _data[id];
    }
    const_reference operator[] (Id const &id) const {
      return (*this)[id.toUnsigned()];
    }

    /// non-const, range checked access into the array
    reference at(Id const &id)
    {
      if (!validId(id)) throw InvalidIndex(id);
      return _data[id.toUnsigned()];
    }

    /// const, range checked access into the array
    const_reference at(Id const &id) const
    {
      if (!validId(id)) throw InvalidIndex(id);
      return _data[id.toUnsigned()];
    }

    /// erase an element -- reset to default value
    void erase (Id const &id)
    {
      if (!validId(id)) throw InvalidIndex(id);
      _data[id.toUnsigned()] = T();
      return;
    }

    /// Delete all the items from the container.
    void clear() {
      tbltype().swap(_data);
    }

    /*!
      compress the array by migrating all values from the old
      location to the new one.  Probably most efficient to build
      new vector then use swap() to trim capacity
    */

    void compress(CompressionData const &cd)
    {
      // stop iterating when we hit the end of the compressed data
      // or the end of this storage.  This allows the sizes of the two
      // objects to mismatch (see CCMP-14338).
      unsigned stopId = std::min(this->maxId().toUnsigned(),
				 cd.getMaxUncompressedId().toUnsigned());

      tbltype newdata(cd.getMaxCompressedId().toUnsigned() + 1);
      
      for (unsigned id = Id::firstValid.toUnsigned(); id <= stopId; ++id) 
	{
	  Id compressedId = cd.getCompressedId(Id(id));
	  if (compressedId.isValid())
	    newdata[compressedId.toUnsigned()] = (*this)[id];
	}
      
      // get rid of excess capacity
      newdata.swap(_data);
    }

  private:
    tbltype _data;
  };       
}


#endif
