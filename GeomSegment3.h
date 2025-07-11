#ifndef _GEOM_SEGMENT3
#define _GEOM_SEGMENT3


#include "GeomBoundedShape3.h"

#include "GeomPoint3.h"

namespace Geom 
{
  /**
   * A geometry class storing a segment in 3-D.  Things like the
   * segment vector and length are cached so subsequent retrieval is
   * faster.
   */
  class Segment3 : public BoundedShape3
  {
  public:
    Segment3()
      :	_hasdirection(false)
      ,	_hasunitdirection(false)
      ,	_haslength(false)
      ,	_hasmidpoint(false)
      , _hasBoundingBox(false)
    {}

    Segment3(Vector3d const (&v)[2])
      : _hasdirection(false)
      ,	_hasunitdirection(false)
      ,	_haslength(false)
      ,	_hasmidpoint(false)
      , _hasBoundingBox(false)
    {
      _vertices[0].position() = v[0];
      _vertices[1].position() = v[1];
    }

    Segment3(Vector3d const &v0,
	     Vector3d const &v1)
      : _hasdirection(false)
      ,	_hasunitdirection(false)
      ,	_haslength(false)
      ,	_hasmidpoint(false)
      , _hasBoundingBox(false)
    {
      _vertices[0].position() = v0;
      _vertices[1].position() = v1;
    }

    Segment3(Segment3 const &rhs)
      : BoundedShape3(rhs)
    {
      *this = rhs;
    }

    Segment3 & operator=(Segment3 const &rhs)
    {
      _vertices[0] = rhs._vertices[0];
      _vertices[1] = rhs._vertices[1];
      if ((_hasdirection = rhs._hasdirection)) _direction = rhs._direction;
      if ((_hasunitdirection = rhs._hasunitdirection)) _unitdirection = rhs._unitdirection;
      if ((_haslength = rhs._haslength)) _length = rhs._length;
      if ((_hasmidpoint = rhs._hasmidpoint)) _midpoint = rhs._midpoint;
      if ((_hasBoundingBox = rhs._hasBoundingBox)) _boundingBox = rhs._boundingBox;
      return *this;
    }


    /*! Get a coordinate of the segment. */
    Vector3d const &getCoordinate(int i) const {
      return _vertices[i].position();
    }

    Vector3d getCoordinate(double p) const {
      Vector3d location;
      location = _vertices[0].position() + (_vertices[1].position()-_vertices[0].position())*p;
      return location;
    }

    /*! Get a point of the segment. */
    Point3 const &getPoint(int i) const {
      return _vertices[i];
    }

    /*! Set the coordinates of the segment. */
    void setCoordinates(Vector3d const &v1,
			Vector3d const &v2) {
      _vertices[0].position() = v1;
      _vertices[1].position() = v2;
      _hasdirection = false;
      _hasunitdirection = false;
      _haslength = false;
      _hasmidpoint = false;
      _hasBoundingBox = false;
    }

    void setIds(Id const &id1,
		Id const &id2) {
      _vertices[0].setId(id1);
      _vertices[1].setId(id2);
    }

    void getIds(Id &id1,
		Id &id2) const {
      id1 = _vertices[0].getId();
      id2 = _vertices[1].getId();
    }
        
    /*! Get the non-unit direction vector (v1->v2) of the segment. */
    Vector3d const &getDirection() const {
      if (!_hasdirection) {
	_direction = _vertices[1].position() - _vertices[0].position();
	_hasdirection = true;
      }
      return _direction;
    }
    /*! Set the non-unit direction vector (v1->v2) of the segment. */
    void setDirection(Vector3d const &direction) {
      _direction = direction;
      _hasdirection = true;
    }

    /*! Get the unit direction vector (v1->v2) of the segment. */
    Vector3d const &getUnitDirection() const {
      if (!_hasunitdirection) {
        this->computeUnitDirection();
      }
      return _unitdirection;
    }
    /*! Set the unit direction vector (v1->v2) of the segment. */
    void setUnitDirection(Vector3d const &unitdirection) {
      _unitdirection = unitdirection;
      // fixme -- add development assertion that vector is unit length
      _hasunitdirection = true;
    }

    /*! Get the length of the segment. */
    double getLength() const {
      if (!_haslength) {
        this->computeUnitDirection(); // need it for length
      }
      return _length;
    }

    /*! Get the length squared -- simple enough, don't bother saving */
    double getLengthSquared() const {
      return this->getCoordinate(0).distanceSquared(this->getCoordinate(1));
    }

    /*! Set the length of the segment. */
    void setLength(double length) {
      _length = length;
      _haslength = true;
    }

    /*! Get the midpoint of the segment. */
    Vector3d const &getMidpoint() const {
      if (!_hasmidpoint) {
	_midpoint = (_vertices[0].position() + _vertices[1].position())/2;
	_hasmidpoint = true;
      }
      return _midpoint;
    }
    /*! Set the midpoint of the segment. */
    void setMidpoint(Vector3d const &midpoint) {
      _midpoint = midpoint;
      _hasmidpoint = true;
    }


    // satisfy BoundedShape3 virtual interface
    BoundingBox getBoundingBox() const;
    BoundingBox const& getConstBoundingBox() const;
    void computeBoundingBox() const;

    bool intersectsBox(BoundingBox const &bbox) const;

    // satisfy Shape3 virtual interface
    Shape3 *clone() const {
      return new Segment3(*this);
    }

    bool containsPoint(Vector3d const &point) const;            

    double getDistanceSquared (Segment3 const &rkSeg1, 
			       double* pfSegP0 = 0, 
			       double* pfSegP1 = 0) const;

    double getDistanceSquared2(Segment3 const &rkSeg1, 
			       double* pfSegP0 = 0, 
			       double* pfSegP1 = 0) const;

    /*! Compute the squared distance from a segment to a point. If the
      argument /a closestpoint is non-null, the closest point on the segment
      will be placed in the input vector. */
    double getDistanceSquared(Vector3d const &querypt,
			      Vector3d *closestpoint = 0,
			      double *param = 0) const;

    /*! Compute the distance from a segment to a point.  If the
      argument /a closestpoint is non-null, the closest point on the segment
      will be placed in the input vector. */
    double getDistance(Vector3d const &point,
		       Vector3d *closestpoint = 0,
		       double *param = 0) const {
      return sqrt(getDistanceSquared(point, closestpoint, param));
    }

  private:
    Point3 _vertices[2];

    mutable bool _hasdirection;
    mutable Vector3d _direction;

    mutable bool _hasunitdirection;
    mutable Vector3d _unitdirection;

    mutable bool _haslength;
    mutable double _length;

    mutable bool _hasmidpoint;
    mutable Vector3d _midpoint;

    void computeUnitDirection() const;

    mutable bool _hasBoundingBox;
    mutable BoundingBox _boundingBox;
  };

}

#endif
