#ifndef _GEOM_STORAGE
#define _GEOM_STORAGE

namespace Geom
{
  class Id;
  class CompressionData;

  /*!
    The Storage class is an abstract base class for data which is accessed by 
    Id.  The base class provides only very rudimentary manipulation
    of the storage so the inefficiencies of polymorphism cannot influence
    critical methods like element access.

    The value of this class is that, along with the StorageGroup class, it
    allows algorithms which modify the mesh to keep mesh-associated data 
    current without explicitly knowing what data is stored.  If data is
    associated with vertices and a vertex is deleted by an algorithm, the
    Storage::erase method can be used to delete the associated data 
    automatically without the algorithm having to know about this data, as 
    long as it is stored in the StorageGroup class.
  */
  class Storage
  {
  public:
    Storage();
    virtual ~Storage() {}

    /*!
      Resize the container to hold Id.  Note, this differs from a normal
      container's resize method in two ways:

      1) It resizes the array to hold the Id passed in, i.e. resize(Id(10)) 
      will actually resize to size 11 so [10] is a valid index.
      2) It will never reduce the size of the array, only grow it.
    */
    virtual void resizeNoShrink(Id const &) = 0;

    /*!
      Resize the container so the maximum Id it can hold is id.   Can shrink
      the array!
    */
    virtual void resize(Id const &) = 0;

    /*!
      Erase the information about Id from storage.
    */
    virtual void erase(Id const &) = 0;

    /*!
      Compress the storage given a CompressionData object.  
    */
    virtual void compress(CompressionData const &cd) = 0;

    /*!
      Clear all data in the container.  Derived classes should implement
      this in a way such that memory is released (i.e. vector<>.swap() instead
      of clear()).  
    */
    virtual void clear() = 0;

    /*!
      Is there data in the storage?
    */
    virtual bool empty() const = 0;
  };
} 

#endif
