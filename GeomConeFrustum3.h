#ifndef _GEOM_CONEFRUSTUM3
#define _GEOM_CONEFRUSTUM3


#include "GeomBoundedShape3.h"
#include "GeomPoint3.h"

#include "GeomBoundingBox.h"

namespace Geom
{
  /*! A ConeFrustum in 3d. */
  class ConeFrustum3 : public BoundedShape3
  {
  public:
    /*! Create the Frustum given the centers and radii of the two circular 
      caps. Internally, the frustum will be stored with the smaller 
      radius end first, so the order of the endpoints may be different than
      that entered. */
    ConeFrustum3(Vector3d const &center1, double radius1,
		 Vector3d const &center2, double radius2);

    /*! Copy constructor. */
    ConeFrustum3(ConeFrustum3 const &rhs);

    /*! Assignment operator. */
    ConeFrustum3 &operator=(ConeFrustum3 const &rhs);

    /*! Clone method. */
    Shape3 *clone() const {
      return new ConeFrustum3(*this);
    }

    /*! Get the center of the Frustum on the side with the smaller radius. */
    Vector3d const &getCenter1() const {
      return _center1;
    }

    /*! Get the smaller radius of the frustum. */
    double const &getRadius1() const {
      return _radius1;
    }

    /*! Get the center of the Frustum on the side with the larger radius. */
    Vector3d const &getCenter2() const {
      return _center2;
    }

    /*! Get the larger radius of the frustum. */
    double const &getRadius2() const {
      return _radius2;
    }

    /*! Implement BoundedShape3::getBoundingBox() method. */
    BoundingBox getBoundingBox() const {
      if (!_hasBoundingBox)
	this->computeBoundingBox();
      return _boundingBox;
    }

    /*! Get the cone angle for the frustum. */
    double getConeAngle() const;

    /*! Get the axis vector (from c2 to c1) for the frustum. 
      (not normalized) */
    Vector3d getAxisVector() const {
      return Vector3d(_center2 - _center1);
    }

    /*! Implement BoundedShape3::intersectsBox() method. */
    bool intersectsBox(BoundingBox const &bbox) const;

    /*! Implement Shape::containsPoint() method. */
    bool containsPoint(Vector3d const &point) const;

    // intersection
    /*! Returns true if the point is contained within the ConeFrustum 
      (inclusive), false otherwise. */
    bool testIntersection(Vector3d const &coord) const;

    bool isDegenerate() const; 

  private:
    Vector3d _center1;
    double _radius1;
    Vector3d _center2;
    double _radius2;
    double _centerLineLength;

    mutable bool _hasBoundingBox;
    mutable BoundingBox _boundingBox;
    void computeBoundingBox() const;

    mutable bool _hasConeVertex;
    mutable Vector3d _coneVertex;
    void computeConeVertex() const;
  };
}
    

#endif
