#ifndef _GEOM_SUBSET_H
#define _GEOM_SUBSET_H

#include <iterator>
#include <ostream>
#include <unordered_set>

#include "GeomId.h"

#include "dynamic_bitset.h"

namespace Geom
{
  // forward declare for friendship
  class Subset;
  class SubsetIterator;
  class CompressionData;

  typedef std::unordered_set<Id, Id::HashFcn, Id::EqualKey> IdHashSet;

  struct SubsetTypes {
    typedef IdHashSet set_type;
    typedef boost::dynamic_bitset < > vector_type;
  };

  /*! A const iterator over subset data.  Iterating over a subset will likely
   *  be faster if the subset is stored with the Set representation.
   *  This iterator should function in any place a std iterator would, except
   *  in modifying algorithms which wouldn't make sense.
   */
  class SubsetIterator :
    public std::iterator<std::forward_iterator_tag, Id, unsigned, Id*, Id &>
  {
  public:
    SubsetIterator() : _subset(0), _setit(), _id() {}
    /*! Construct the iterator given a const pointer to a subset. */
    explicit SubsetIterator(Subset const *subset)
      : _subset(subset), _setit(), _id() {
      this->begin();
    }
    // use default copy constructor and assignment

    void begin();			// set to first id in set
    bool done() const			{ return !_id.isValid(); }
    SubsetIterator &operator++();

    Id const& operator*() const		{ return _id; }
    Id const *operator->() const	{ return &_id; }

    bool operator!=(SubsetIterator const &si)  const { return _id != si._id; }
    bool operator==(SubsetIterator const &si)  const { return _id == si._id; }

  protected:
    Subset const *_subset;
    SubsetTypes::set_type::const_iterator _setit;
    Id _id;				// current id, invalid if at end

    // since most people use postfix increment without realizing the added cost,
    // it is best to block its usage, and there isn't a strong need for it.
    SubsetIterator operator++(int);
  };

  /*! The Subset class is used to store topological subsets such as
    cell sets, face sets (surfaces), etc.  Internally, it allows these
    subsets to be stored using one of two different representations, a 
    set representation and a vector representation.  Each of these
    representations has advantages and disadvantages, so the proper
    representation should be chosen for each algorithm.

    The Set representation uses an associative container (hash_set)
    to store the elements in the subset.  This representation is best
    for algorithms for which the subset is likely to be small compared
    to the entire set (i.e. boundary faces are small compared to the number
    of interior faces).  This representation is also best when it is more 
    important to quickly iterate over the elements of the subset rather than 
    quickly determine if an individual element is in the subset.
        
    The Vector representation uses a vector of boolean values
    to store the elements of the subset.  This representation is best
    for algorithms for which the subset is not small compared to the
    size of the set (i.e. a cell set which could contain all cells).  
    This representation provides fast determination of whether an individual
    element is in the subset, but iteration over the elements of the
    subset is slower than the set representation.
  */
         
  class Subset 
  {
  public:
    enum Mode {
      Set = 0x01,
      Vector = 0x02
    };

    typedef Id value_type;
    typedef SubsetTypes::set_type set_type;
    typedef SubsetTypes::vector_type vector_type;
    typedef SubsetIterator iterator;
    typedef SubsetIterator const_iterator;

    char const * modeToChar(Mode mode)
      { return (mode == Set) ? "Set" : "Vector"; }

    /*! Construct the subset.  By default, the set representation is 
      chosen initially. */
    Subset(Mode mode = Set);

    template<typename FORWARD_ITERATOR>
    Subset(FORWARD_ITERATOR const & beg,
                   FORWARD_ITERATOR const & end, Mode mode = Set)
      : _mode(mode)
      , _minid(Id::max)
      , _maxid(Id::min)
      , _numelements(0)
    { add(beg, end); }

    template<typename CONTAINER>
    explicit Subset(CONTAINER const & cont, Mode mode = Set)
      : _mode(mode)
      , _minid(Id::max)
      , _maxid(Id::min)
      , _numelements(0)
    { add(cont); }

    // copy constructor and assignment operator
    Subset(Subset const &rhs);
    Subset &operator=(Subset const &);

    ~Subset();
        
    /*! Set the representation mode of the subset.  If the mode is the
      same as the current mode, nothing will happen.  Otherwise, the
      subset data will be converted into the new representation and 
      the subset data from the previous representation will be cleared.
    */
    void setMode(Mode mode);
    /*! Retrieve the representation mode of the subset. */
    Mode getMode() const {return _mode;}
        
    /*! Clear the subset.  This operation will not resize any data but
      will remove all elements from the subset. */
    void clear();

    /*! Really clear the subset.  This operation will also free all storage
      associated with the subset. */
    void clearAll();

    /*! Resize the data to accomatate the specified maximum id.  The
      data will grow or shrink to accomodate added Id's, but it is better
      to explicitly resize the data whenever possible. */
    void resize(Id const &maxid);

    /*! Add an element to the subset. If the element is already in the
      subset, nothing will happen. */
    bool add(Id const &id);

    /*! Add a subset to the subset.  This is equivalent taking the union
      of two subsets.
    */
    void add(Subset const &subset);

    /*! Add all ids from iterator range */
    template<typename FORWARD_ITERATOR>
    void add(FORWARD_ITERATOR const &beg, FORWARD_ITERATOR const &end) {
      for (FORWARD_ITERATOR i(beg); i != end; ++i)
	this->add(*i);
    }

    /*! Add all ids from any container that defines begin, end,
     *  and a forward iterator.
     */
    template<typename CONTAINER>
    void add(CONTAINER const &cont) {
      add(cont.begin(), cont.end());
    }

    /*! Intersect two subsets.  After this call, all elements that are in
      this subset AND in the input subset will be in this subset.   Any 
      elements in this but not in the input subset will be removed.
    */
    void intersect(Subset const &subset);

    /*! Remove an element from the subset.  If the element is not contained
      in the subset, nothing will happen.  Returns true if an element was 
      removed from the subset. */
    bool remove(Id const &id);

    /*! Remove a subset from the subset. Returns true if an element was 
      removed from the subset. */
    bool remove(Subset const &subset);

    /*! Remove all ids from iterator range from the subset.
      Returns true if an element was removed from the subset. */
    template<typename FORWARD_ITERATOR>
    bool remove(FORWARD_ITERATOR const &beg, FORWARD_ITERATOR const &end) {
      bool anyRemoved = false;
      for (FORWARD_ITERATOR i(beg); i != end; ++i)
	if (this->remove(*i))
          anyRemoved = true;
      return anyRemoved;
    }

    /*! Remove container of ids from the subset. 
      Returns true if an element was removed from the subset. */
    template<typename CONTAINER>
    bool remove(CONTAINER const &cont) {
      return remove(cont.begin(), cont.end());
    }

    /*! Check if the subset contains a given element. This method will
      work faster if the subset is stored with the Vector representation. */
    bool contains(Id const &id) const;              

    /*! Return the minimum Id stored in the subset */
    Id getMinId() const;

    /*! Return the maximum Id stored in the subset */
    Id getMaxId() const;

    /*! Return true if the subset is empty, else false. */
    bool empty() const;

    /*! Return the number of elements stored in the subset. */
    unsigned getNumElements() const;

    // begin/end iterators
    SubsetIterator begin() const        { return SubsetIterator(this); }
    SubsetIterator end() const          { return SubsetIterator(); }

    void swap(Subset &subset);

    /*! Compress the subset, generating a new instance of 
      CompressionData.  After completion, the Subset will have
      the following properties:

      * getMinId() == Id::firstValid
      * getMaxId() == getNumElements
      * Be of dense (vector) type, since it is by default dense.

      The mapping from compressed to uncompressed Ids will be stored in the
      CompressionData object. This mapping must then be passed into all
      other subsets in the same Id space to keep Ids synchronized.  */
    void compressMaster(CompressionData &data);

    /*! Compress the subset based on an existing CompressionData.  
      Call this method when another Subset referring to the same Id space
      has been compressed using Subset::compressMaster.

      Note that it only makes sense for all slaves to be proper subsets of the
      master subset for which compressMaster was called.  If any Ids are found
      in the slave subset that were not compressed in the master subset, 
      Id collisions are possible.  Thus, an exception will be thrown if an
      Id is found in the slave that was not in the master.

      This note along with the assumption that all current subsets on a 
      particular Id space will be compressed implies that the subset for which
      compressMaster is called must be exactly the union of all other subsets
      on that id space.
    */
    void compressSlave(CompressionData const &data);

    /**
     * Compare this subset with another for strict equality.  All elements 
     * in this subset must be in the other one and vice versa.
     */
    bool equals(Subset const &rhs) const;

  private:
    Mode _mode;

    mutable Id _minid;
    mutable Id _maxid;
    mutable unsigned _numelements;

    set_type _setdata;
    vector_type _vectordata;

    friend class SubsetIterator;

    friend std::ostream & operator<<(std::ostream &s, Subset const &subset);
  };

  /** Insert iterator for Subsets.
   *  Can be used to insert items from another container into a Subset.
   *  For example, to copy the contents of some container into a mySubset do:
   *  copy(cont.begin(), cont.end(), SubsetInsertIterator(mySubset));
   */
  class SubsetInsertIterator
    : public std::iterator<std::output_iterator_tag, void, void, void, void>
  {
  private:
    Subset * const _subset;
  public:
    typedef Subset container_type;
    explicit SubsetInsertIterator(Subset &ss) : _subset(&ss) {}
    SubsetInsertIterator & operator=(Id id) {
      _subset->add(id);
      return *this;
    }
    SubsetInsertIterator & operator=(SubsetInsertIterator const & rhs)
    {
      if (this != &rhs)
        const_cast<Subset *&>(_subset) = rhs._subset;
      return *this;
    }
    SubsetInsertIterator & operator*()     { return *this; }
    SubsetInsertIterator & operator++()    { return *this; }
    SubsetInsertIterator & operator++(int) { return *this; }
  };

}

#endif
