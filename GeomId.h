#ifndef _GEOM_ID_H
#define _GEOM_ID_H

#include <iosfwd>

#include <type_traits>

namespace Geom
{
  class Id
  {
  public:
    // use these typedefs instead of int when converting an Id to an int type.
    typedef int SignedType;
    typedef unsigned int UnsignedType;

    //! min is an Id such that min < id for all valid ids
    static Id const min;
    //! max is an Id such that max > id for all valid ids
    static Id const max;
    //! firstValid is the minimum id > min
    static Id const firstValid;
    
    Id()
      : _id(max._id)
    {}

    explicit Id(UnsignedType id, bool sign = true)
      : _id(sign ? id : -id)
    {}

    void setId(Id id) { _id = id._id; }
    Id set(UnsignedType value, bool sign = true)
    {
      _id = sign ? value : -value;
      return *this;
    }

    bool getSign() const {
      return (_id >= 0);
    }

    void setSign(bool sign) {
      if (sign != getSign())
        _id = -_id;
    }

    void invertSign() {
        _id = -_id;
    }

    UnsignedType getValue() const {
      return (getSign() ? _id : -_id);
    }

    void setValue(UnsignedType id) {
      _id = (getSign() ? SignedType(id) : -SignedType(id));
    }

    UnsignedType toUnsigned() const {
      return getValue();
    }

    bool isValid() const {
      return (_id != 0 && getValue() < max.getValue());
    }

    /*! Added this operator so we can do !id instead of !id.isValid(). */
    bool operator!() const {
      return !isValid();
    }

    void swap(Id &x) {
      SignedType const xid = x._id;
      x._id = _id;
      _id = xid;
    }

    bool operator==(Id const &b) const {
      return (getValue() == b.getValue());
    }

    bool identical(Id const &b) const
      { return (_id == b._id); }
  
    bool validAndIdentical(Id const &b) const
      { return (isValid() && b.isValid() && _id == b._id); }

    bool operator!=(Id const &b) const {
      return (getValue() != b.getValue());
    }
  
    bool operator<(Id const &b) const {
      return (getValue() < b.getValue());
    }
  
    bool operator<=(Id const &b) const {
      return (getValue() <= b.getValue());
    }
  
    bool operator>(Id const &b) const {
      return (getValue() > b.getValue());
    }
  
    bool operator>=(Id const &b) const {
      return (getValue() >= b.getValue());
    }

    /*! A Hash Functor for the Id.  Ids are hashed with a comparison function
      that ignores their sign. */
    struct HashFcn
    {
      typedef UnsignedType return_type;
      inline return_type operator() (Id const &key) const {
        return key.getValue();
      }
    };
    
    /*! A comparison functor for the Id.  Ids are compared ignoring their 
      sign. */
    struct EqualKey
    {
      inline bool operator() (Id const &a, Id const &b) const {
        return (a.getValue() == b.getValue());
      }
    };


    /*! A comparison functor for the Id.  Ids are compared including their 
      sign. */
    struct IdenticalKey
    {
      inline bool operator() (Id const &a, Id const &b) const {
	return (a._id == b._id);
      }
    };

    struct IdenticalStrictWeakOrdering
    {
      bool operator()(Id const &a, Id const &b) const
      { return (a._id < b._id); }
    };


    /**
     * Return Id that is different, behavior undefined if all different.
     * For example, if a == c, then getUnique(a, b, c) == b.
     *
     * The getUnique functions are useful for getting the Id that is 
     * different.  For example, if you have the three vertices on a
     * triangular face (v0,v1,v2), and the two vertices on one of its 
     * edges (v3,v4), then getUnique(v0,v1,v2,v3,v4) will return the
     * vertex not on the edge.
     * These are implemented using exclusive or (xor) as this doesn't
     * suffer from overflow problems that using addition has.
     */
    static Id getUnique(Id a, Id b, Id c) {
      return Id(a.getValue() ^ b.getValue() ^ c.getValue());
    }

    /**
     * Return Id that is different, behavior undefined if more than one unique.
     * For example, if a == c and b == e, then getUnique(a, b, c, d, e) == d.
     */
    static Id getUnique(Id a, Id b, Id c,
			Id d, Id e) {
      return Id(a.getValue() ^ b.getValue() ^ c.getValue() ^ d.getValue() ^ e.getValue());
    }


  private:
    SignedType _id;
  };


  std::ostream &operator<< (std::ostream &s, Id const &id);
  std::istream &operator>> (std::istream &s, Id       &id);

  typedef Id VertId;
  typedef Id EdgeId;
  typedef Id FaceId;
  typedef Id CellId;
}

#endif
