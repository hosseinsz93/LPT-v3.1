#ifndef _GEOM_BOX_H
#define _GEOM_BOX_H

#include <cmath>
#include <iosfwd>

#include "GeomVector.h"

namespace Geom
{
  /*! A BoundingBox class in 3d. */
  class BoundingBox
  {
  public:
    enum { SIZE = 2 };          // number of corners

    /*! Construct an empty box.  The box will have min values set
      to a large positive value and max values set to a large minimum
      value so the first call to "expand" will initialize it. */
    BoundingBox() { reset(); }

    /*! Construct an uninitialized box.*/
    BoundingBox(UninitializedToken const &) {}

    // Use default copy constructor and assignment operator.

    /*! Construct a box around a point. */
    explicit BoundingBox(Vector3d const &coord) {
      reset(coord);
    }

    /*! Construct a box with specified min and max points. */
    BoundingBox(Vector3d const &minpnt,
                Vector3d const &maxpnt) {
      reset(minpnt, maxpnt);
    }

    /*! Return the lower corner (smallest x, y, z). */
    Vector3d const &min() const { return _min(); }
    /*! Return the upper corner (largest x, y, z). */
    Vector3d const &max() const { return _max(); }

    /*! Return corner n (must be zero or one). */
    Vector3d const &corner(int n) const {
      return _corner[n];
    }
    /*! Return (i,j,k) corner (i,j,k must be zero or one). */
    Vector3d corner(int i, int j, int k) const {
      return Vector3d(_corner[i](0), _corner[j](1), _corner[k](2));
    }

    /*! Reset the bounding box to default settings (large min, small max) */
    void reset();
    void reset(Vector3d const &point);
    void reset(Vector3d const &minc, Vector3d const &maxc);

    /*! Expand the box to hold the specified point. */
    void expand(Vector3d const &coord);

    /*! Expand the box to hold the box. */
    void expand(BoundingBox const &bb);

    /*! Expand the box to hold points in container range. */
    template<typename ITERATOR>
      void expand(ITERATOR const& begin, ITERATOR const& end) {
        for (ITERATOR i(begin); i != end; ++i)
          expand(*i);
    }

    /*! Inflate the box by the specified amount in all directions.
     *  The total inflation amount is 2*size in each direction */
    void inflate(double const size);

    /*! Inflate the box by the specified amount in each directions.
     *  The total inflation amount is 2*size[] in each direction */
    void inflate(Vector3d const &size);

    /*! Check if the specified point is contained in the box (inclusive) */
    bool containsPoint(Vector3d const &coord) const {
      // optimizeme:  try other checks here...
      return (coord(0) >= _min(0) &&
              coord(1) >= _min(1) &&
              coord(2) >= _min(2) &&
              coord(0) <= _max(0) &&
              coord(1) <= _max(1) &&
              coord(2) <= _max(2));
    }

    /*! Check if this box intersects with the specified box. */
    bool intersects(BoundingBox const &rhs) const;

    /*! Check if this box contains the specified box. */
    bool contains(BoundingBox const &rhs) const;

    /*! Check if this box intersects with the specified triangle . */
    bool intersectsTriangle(Vector3d const &v1,
                            Vector3d const &v2,
                            Vector3d const &v3,
                            double tolerance = 1.0e-9) const;

    /*! Check if this box intersects with the specified sphere. */
    bool intersectsSphere(Vector3d const &point,
                          double radius) const;

    /*! Check if this box intersects with the given ray. */
    bool intersectsRay(Vector3d const &point,
                       Vector3d const &dir,
                       double *distance) const;

    /*! use intersectsRay to determine if segment (p1,p2) intersects box */
    bool intersectsSegment(Geom::Vector<3, double> const &p1,
                           Geom::Vector<3, double> const &p2);

    /*! Compute the minimum distance between the two boxes. */
    double distance(BoundingBox const &bb) const {
      return sqrt(distanceSquared(bb));
    }

    /*! Compute the minimum distance between the box and a point.
      If the point lies inside the box, returns zero. */
    double distance(Vector3d const &coord) const {
      return sqrt(distanceSquared(coord));
    }

    /*! Compute the minimum squared distance between the two boxes. */
    double distanceSquared(BoundingBox const &bb) const;

    /*! Compute the minimum squared distance between the box and a point.
      If the point lies inside the box, returns zero. */
    double distanceSquared(Vector3d const &coord) const;

    /*! Return center of the box. */
    Vector3d center() const {
      return Vector3d(0.5*(_min(0)+_max(0)),
                      0.5*(_min(1)+_max(1)),
                      0.5*(_min(2)+_max(2)));
    }

    /*! Return major diagonal. */
    Vector3d diagonal() const {
      return Vector3d(sideLength(0), sideLength(1), sideLength(2));
    }

    /*! Return boxes ith dimension. */
    double sideLength(int i) const {
      return _max(i) - _min(i);
    }

    /*! Return the longest dimension of the box. */
    double longestSideLength() const {
      return std::max(std::max(sideLength(0), sideLength(1)), sideLength(2));
    }

  protected:

    // low-level accessors
    Vector3d const &_min() const { return _corner[0]; }
    Vector3d const &_max() const { return _corner[1]; }
    Vector3d &_min()             { return _corner[0]; }
    Vector3d &_max()             { return _corner[1]; }
    double const &_min(int i) const      { return _corner[0](i); }
    double const &_max(int i) const      { return _corner[1](i); }
    double &_min(int i)                  { return _corner[0](i); }
    double &_max(int i)                  { return _corner[1](i); }

  private:

    Vector3d _corner[SIZE];     // min/max corners of box

  };

  std::ostream & operator<< (std::ostream &out, BoundingBox const &bbox);

}

#endif
