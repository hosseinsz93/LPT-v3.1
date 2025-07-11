#include "GeomSubset.h"
#include "GeomCompressionData.h"
#include "GeomUtility.h"

#include <algorithm>

namespace Geom
{
  Subset::Subset(Mode mode)
    : _mode(mode),
      _minid(Id::max),
      _maxid(Id::min),
      _numelements(0)
  { }
    
  Subset::~Subset()
  { }

  Subset::Subset(Subset const &rhs)
    : _mode(rhs._mode), 
      _minid(rhs._minid),
      _maxid(rhs._maxid),
      _numelements(rhs._numelements),
      _setdata(rhs._setdata),         
      _vectordata(rhs._vectordata)
  { }

  Subset &Subset::operator=(Subset const &rhs)
  {
    if (this == &rhs)
      return *this;

    _mode = rhs._mode;
    _maxid = rhs._maxid;
    _minid = rhs._minid;
    _numelements = rhs._numelements;
    _setdata = rhs._setdata;
    _vectordata = rhs._vectordata;
    return *this;
  }

  void Subset::setMode(Mode mode)
  {
    if (_mode == mode)
      return; // in same mode

    if (mode == Set)
      {
        // switching to set mode
        _setdata.clear();

        vector_type::size_type pos;
        for (pos=_vectordata.find_first();
	     pos!=_vectordata.npos();
	     pos=_vectordata.find_next(pos))
          _setdata.insert(Id(pos));

        // clear vector mode data using "swap trick"
        vector_type().swap(_vectordata);
      }
    else
      {
        // switching to vector mode
        Id maxid = this->getMaxId();
        if (maxid != Id::min)
          _vectordata.resize(maxid.toUnsigned() + 1, false);
        _vectordata.reset();
        for (set_type::const_iterator it = _setdata.begin(); 
	    it != _setdata.end(); ++it)
          _vectordata[it->toUnsigned()] = true;
                            
        // clear set data
        _setdata.clear();
      }


    _mode = mode;   
  }

  void Subset::clear()
  {
    if (_mode == Set)
      _setdata.clear();
    else
      _vectordata.reset();

    _minid = Id::max;
    _maxid = Id::min;
    _numelements = 0;
  }

  void Subset::clearAll()
  {
    _setdata.clear();
    vector_type().swap(_vectordata);

    _minid = Id::max;
    _maxid = Id::min;
    _numelements = 0;
  }

  void Subset::resize(Id const &maxid)
  {
    if (_mode == Set)
      return;     // nothing to do in set mode

    if (maxid.toUnsigned() >= _vectordata.size())
      _vectordata.resize(maxid.toUnsigned() + 1, false);
  }
            

  bool Subset::add(Id const &id)
  {
    if (_mode == Set)
      {
	// Note: Id::HashFcn converts id to positive so sign is ignored.
        if (!_setdata.insert(id).second)
          return false;		// already in set

	if (_numelements == 0)
	  {
	    // first element in set -- obviously min and max id
	    _numelements = 1;
	    _minid = id;
	    _maxid = id;
	  }
	else
	  {
	    ++_numelements;

	    // only update min and max id if they're currently valid -- 
	    // otherwise, save for later.
	    if (_minid != Id::max)
	      _minid = std::min(_minid, id);
	    if (_maxid != Id::min)
	      _maxid = std::max(_maxid, id);
	  }
      }
    else
      {
        this->resize(id);
        if (_vectordata[id.toUnsigned()])
          return false;		// already in set
        _vectordata[id.toUnsigned()] = true;

	++_numelements;
	// update min and max id, since this can be done quickly
	_minid = std::min(_minid, id);
	_maxid = std::max(_maxid, id);
      }

    return true;
  }

  void Subset::add(Subset const &subset)
  {
    if (&subset == this)
      return; // nothing to do

    // for now, just use the subset iterator and existing add method.
    // this could be made more efficient at some point
    for (SubsetIterator si(&subset); !si.done(); ++si)
      this->add(*si);
  }

  void Subset::intersect(Subset const &subset)
  {
    if (&subset == this)
      return; // nothing to do

    Subset elementsToRemove;
    for (SubsetIterator si(this); !si.done(); ++si)
      if (!subset.contains(*si))
	elementsToRemove.add(*si);
    remove(elementsToRemove);
  }
  
  bool Subset::remove(Id const &id)
  {   
    if (_mode == Set)
      {
	set_type::iterator it = _setdata.find(id);
	if (it != _setdata.end())
	  {
	    _setdata.erase(it);
	    --_numelements;
	    
	    // reset min, max to force recomputation later
	    _minid = Id::max;
	    _maxid = Id::min;
	    return true;
	  }
	else
	  return false;
      }
    else
      {
        if (id.toUnsigned() >= _vectordata.size())
          return false;
	if (!_vectordata[id.toUnsigned()])
	  return false;

        _vectordata[id.toUnsigned()] = false;

        // automatically keep track of min and max id in vector
        // rep since these are relatively cheap operations
        _numelements--;
        if (_numelements == 0)
          {
            _minid = Id::max;
            _maxid = Id::min;
          }
        else
          {
            if (id == _minid)
              {
                for (unsigned i = _minid.toUnsigned(); 
                    i < _vectordata.size(); ++i)
                  if (_vectordata[i])
                    {
                      _minid = Id(i);
                      break;
                    }
              }
            else if (id == _maxid)
              {
                for (unsigned i = _maxid.toUnsigned();
                    i > 0; i--)
                  if (_vectordata[i-1])
                    {
                      _maxid = Id(i-1);
                      break;
                    }
              }
          }
      }

    return true;
  }   


  bool Subset::remove(Subset const &subset)
  {
    if (&subset == this)
      {
	this->clear();
	return true;
      }

    //there are times when this subset can be much smaller than
    //the subset that is passed in.  If we iterate over the smaller
    //subset it can improve performance.  There may be some things
    //that could be done at the stl level to further improve performance
    //but this helps quite a bit.
    if (this->getNumElements() >= subset.getNumElements())
    {
    // for now, just use the subset iterator and existing remove method.
    // this could be made more efficient at some point
      bool removed = false;
      for (SubsetIterator si(&subset); !si.done(); ++si)
        if (this->remove(*si))
          removed = true;
      return removed;
    }
    else
    {
      bool removed = false;
      for (SubsetIterator si(this); !si.done(); )
      {
        Id element(*si);
        ++si;           //increment before removing so interator stays valid
        if (subset.contains(element))
        {
          removed = true;
          this->remove(element);
        }
      }
      return removed;
    }
  }    

    
  bool Subset::contains(Id const &id) const
  {
    if (_mode == Set)
      return (_setdata.find(id) != _setdata.end());
    else
      {
        if (id.toUnsigned() >= _vectordata.size())
          return false;
        return _vectordata[id.toUnsigned()];
      }
  }

  Id Subset::getMinId() const
  {
    if (_minid != Id::max)
      return _minid;

    if (_mode == Set)
      {
        set_type::const_iterator it = std::min_element(_setdata.begin(), 
							_setdata.end());
        if (it != _setdata.end())
          _minid = (*it);
      }
    else
      {
        vector_type::size_type pos = _vectordata.find_first();
        if (pos != vector_type::npos())
          _minid = Id(pos);
      }
    return _minid;
  }

  Id Subset::getMaxId() const
  {
    if (_maxid != Id::min)
      return _maxid;

    if (_mode == Set)
      {
        set_type::const_iterator it = std::max_element(_setdata.begin(), 
							_setdata.end());
        if (it != _setdata.end())
          _maxid = (*it);
      }
    else
     for (unsigned i = _vectordata.size(); i > 0; --i)
        if (_vectordata[i-1])
          {
            _maxid = Id(i-1);
            break;
          }
    return _maxid;
  }

  bool Subset::empty() const
  {
    if (_mode == Set)
      return _setdata.empty();
    // else Vector - _numelements is kept up to date
    return _numelements == 0;
  }

  unsigned Subset::getNumElements() const
  {
    if (_numelements != 0)
      return _numelements;

    if (_mode == Set)
      _numelements = _setdata.size();
    else
      _numelements = _vectordata.count();

    return _numelements;
  }

  void Subset::swap(Subset &subset)
  {
    std::swap(_mode, subset._mode);
    std::swap(_minid, subset._minid);
    std::swap(_maxid, subset._maxid);
    std::swap(_numelements, subset._numelements);
    _setdata.swap(subset._setdata);
    _vectordata.swap(subset._vectordata);
  }

  void Subset::compressMaster(CompressionData &data)
  {
    // cache current number of elements and maximum uncompressed id
    Id numElements(this->getNumElements());
    Id maxId(this->getMaxId());    

    // initialize data
    data.initialize(maxId, numElements);

    // build compressed->uncompressed mapping
    unsigned currentId = Id::min.toUnsigned();
    for (SubsetIterator si(this); !si.done(); ++si) {
      ++currentId;
      data.setCompressedId(*si, Id(currentId));
    }
    
    // clear set data if needed
    if (_mode == Set) {
      _setdata.clear();
      _mode = Vector;
    }

    // resize vectordata to new size.  Use swap instead of resize to
    // reduce capacity
    vector_type(numElements.toUnsigned()+1).swap(_vectordata);
    _vectordata.set();
    _vectordata[0] = false;
    
    // update min and max id
    _maxid = numElements;
    _minid = Id(1);
  }
  
  void Subset::compressSlave(CompressionData const &data)
  {
    // initialize min and max ids
    _minid = Id::max;
    _maxid = Id::min;
    
    // handle set rep and vector rep separately 
    
    if (_mode == Set)
      {
	// create new set, populate it, and swap it with existing set
	set_type newset;

	for (set_type::iterator it = _setdata.begin(); it != _setdata.end();
	    ++it)
	  {
	    Id compressedId = data.getCompressedId(*it);
	    _minid = std::min(_minid, compressedId);
	    _maxid = std::max(_maxid, compressedId);

	    newset.insert(compressedId);
	  }
	  
	
	// swap old and new vectors
	newset.swap(_setdata);
      }
    else
      {
	// for now, create new vector, populate it, and swap it with existing
	// vector
	vector_type newvector(data.getMaxCompressedId().toUnsigned()+1);
        newvector.reset();

	// use subset iterator since it takes care of operator++
	for (SubsetIterator si(this); !si.done(); ++si)
	  {
	    Id compressedId = data.getCompressedId(*si);
	    _minid = std::min(_minid, compressedId);
	    _maxid = std::max(_maxid, compressedId);

	    newvector[compressedId.toUnsigned()] = true;
	  }

	// swap old and new vectors
	newvector.swap(_vectordata);	
      }
  }


  bool Subset::equals(Subset const &rhs) const
  {
    if (rhs.getNumElements() != getNumElements())
      return false;

    for(SubsetIterator si(this); !si.done(); ++si)
      {
	if (!rhs.contains(*si))
	  return false;
      }
    return true;
  }

  std::ostream &operator<<(std::ostream &s, Subset const &subset)
  {
    SubsetIterator si(subset.begin());
    if (!si.done())
      {
	s << (*si);
	for (++si; !si.done(); ++si)
	  s << ',' << (*si);
      }
    return s;
  }


  void SubsetIterator::begin()
  {
    if (!_subset)
      {
	_id = Id::max;
      }
    else if (_subset->_mode == Subset::Set)
      {
	_setit = _subset->_setdata.begin();
	_id = (_setit != _subset->_setdata.end()) ? *_setit : Id::max;
      }
    else
      {
	Subset::vector_type::size_type pos = _subset->_vectordata.find_first();
	_id = (pos != _subset->_vectordata.npos()) ? Id(pos) : Id::max;
      }
  }


  SubsetIterator &SubsetIterator::operator++()
  {
    if (_subset->_mode == Subset::Set)
      {
	_setit++;
	_id = (_setit != _subset->_setdata.end()) ? *_setit : Id::max;
      }
    else
      {
	Subset::vector_type::size_type pos = _id.toUnsigned();
	pos = _subset->_vectordata.find_next(pos);
	_id = (pos != _subset->_vectordata.npos()) ? Id(pos) : Id::max;
      }
    return *this;
  }


  SubsetIterator SubsetIterator::operator++(int)
  {
    SubsetIterator save(*this);
    ++(*this);
    return save;
  }
        
}
