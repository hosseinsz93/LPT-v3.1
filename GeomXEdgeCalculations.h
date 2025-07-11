#ifndef _GEOMX_EDGECALCULATIONS_H_
#define _GEOMX_EDGECALCULATIONS_H_


#include <vector>

#include "GeomVector.h"

#include "GeomEdgeCalculations.h"
// #include "GeomMeshContainer3d.h"
// #include "GeomXVertexCoordInterface.h"

#include "GeomXCross_t.h"


namespace GeomX
{
  
  // get the endpoint coordinates of an edge
#if 0
  void edgeCoordinates(
    MKXCore::VertexCoordInterface * geometry, MKCore::TopologyData *topology,
    MKCore::Id const &edge, Vector3d &coord1, Vector3d &coord2);
#endif
  
  bool edgeEdgeIntersect(Vector2d const &p0,
                         Vector2d const &p1,
                         Vector2d const &q0,
                         Vector2d const &q1,
                         double &alphap, double &alphaq);


  bool edgeEdgeIntersect(Vector2d p0, Vector2d p1,
                         Vector2d q0, Vector2d q1,
                         double linearTol = 1.0e-14);


  Vector3d edgeEdgeIntersect(Vector3d p0, Vector3d p1,
                             Vector3d q0, Vector3d q1);


  template < typename T >
  double findInterpFactor(T a, T p0, T p1)
  {
    T v0 = T(p1 - p0);
    T v1 = T(a  - p0);
    return v0.dot(v1) / v0.mag2();
  }

  template < >
  inline double findInterpFactor(double a, double p0, double p1)
  {
    double v0 = p1 - p0;
    double v1 = a  - p0;
    return v1 / v0;
  }


  template < typename T >
  T interpEdgePosition(unsigned seg0, double coef,
                       std::vector < T > const & coords)
  {
    unsigned seg1 = seg0 + 1;
    if (seg1 >= coords.size())
      seg1 = 0;
    T result(coords[seg0] * (1.0-coef) + coords[seg1] * coef);
    return result;
  }



  // Routines after this comment were added around 2009-Feb-25

  // Ignore boundaries of the segment.
  template < int N, typename T >
  double linePointDistance(VectorNT const &coord1,
                           VectorNT const &coord2,
                           VectorNT const &querypt);

  // Find the closest point on a line to the given point.
  // p0 is a point on the line.  pv is the direction of the line.
  // The return value, ratio, has the value such the p0+ratio*pv is the
  // closest point.

  // Line is used in the title since they are assumed to be infinitely long
  // and the closest point doesn't have to be on the line segment.
  template < int N, typename T > T
  lineProjectPointVec_t(VectorNT const & pnt,
                        VectorNT const & p0, VectorNT const & pv);

  template < int N, typename T >
  inline T lineProjectPointPnt_t(VectorNT const & pnt,
                                 VectorNT const & p0, VectorNT const &p1)
    { return lineProjectPointVec_t(pnt, p0, VectorNT(p1 - p0)); }


  // Line is used in the title since they are assumed to be infinitely long
  // and the intersection point doesn't have to be on the line segment.  If
  // the lines are parallel, only one point is returned which will be located
  // at one edge of one of the segments.
  template < typename T > void
  lineLineIntersectPnt_t(Geom::Vector<2,T> const & p0, Geom::Vector<2,T> const & p1,
                         Geom::Vector<2,T> const & q0, Geom::Vector<2,T> const & q1,
                         T & ratiop, T & ratioq);


  template < typename T > void
  lineLineIntersectPnt_t(Geom::Vector<3,T> const & p0, Geom::Vector<3,T> const & p1,
                         Geom::Vector<3,T> const & q0, Geom::Vector<3,T> const & q1,
                         T & ratiop, T & ratioq);


  // Find the distance squared between an edge (line segment) and a point.
  // This calculation takes into account the length of the segment so that the
  // minimum distance is the distance between pnt and the end points p0 and p1.
  // If the length of p1-p0 is < tol, the edge is assumed to be a point.  It
  // also can calculate the cross product between p1-p0 cross pnt-p0 which can
  // be used to see if points are on the same or different sides of the
  // edge p1-p0.
  template < int N, typename T >
  T edgePointDistance2Pnt_t(VectorNT const & pnt,
                            VectorNT const & p0, VectorNT const & p1,
                            T tol, Cross_t<N,T> * cross = 0);

  template < int N, typename T >
  inline T edgePointDistancePnt_t(VectorNT const & pnt,
                                  VectorNT const &p0, VectorNT const &p1,
                                  T tol = 1.0e-15)
    { return std::sqrt(edgePointDistance2Pnt_t(pnt, p0, p1, tol)); }


  // Find the vector from p0 and p1 to pnt using the edge as the x axis of
  // a coordinate system.  The other axes don't matter since the other result
  // is the normal distance to the line.  The first coordinate in the resulting
  // vector is the distance along the edge from p[01] and the second coordinate
  // is the distance from the line.
  template < int N, typename T >
  void edgePointVectors(VectorNT const & pnt,
                        VectorNT const & p0,     VectorNT const & p1,
                        VectorNT       & fromp0, VectorNT       & fromp1);


  // Find the intersection point of two edges.  In 3d, this would be the
  // closest distance.  If the line segments don't intersect, the routine
  // returns the closest points between the segments.  The points are
  // returned as ratios along the vector from point 0 to point 1.
  template < int N, typename T >
  void edgeEdgeIntersectPnt_t(VectorNT const &p0, VectorNT const &p1,
                              VectorNT const &q0, VectorNT const &q1,
                              T & ratiop, T & ratioq);


  template < int N, typename T >
  T edgeEdgeDistancePnt_t(VectorNT const & p0, VectorNT const & p1,
                          VectorNT const & q0, VectorNT const & q1,
                          T tol);




  template < int N, typename T >
  inline VectorNT findProjectedPointPnt_t(T const & ratio,
                                      VectorNT const & p0, VectorNT const & p1)
  {
    VectorNT pv(p1 - p0);
    return findProjectedPointVec_t(ratio, p0, pv);
  }

  template < int N, typename T >
  inline VectorNT findProjectedPointVec_t(T const & ratio,
                                      VectorNT const & p0, VectorNT const & pv)
  {
    return VectorNT(p0 + ratio * pv);
  }

  template < int N, typename T >
  inline T findProjectedDistancePnt_t(VectorNT const & pnt, T const & ratio,
                                      VectorNT const & p0, VectorNT const & p1)
  {
    VectorNT projPnt(findProjectedPointPnt_t(ratio, p0, p1));
    return VectorNT(projPnt - pnt).length();
  }




  template < int N, typename T >
  inline VectorNT interpolate(double coef, VectorNT const & p0,
                                           VectorNT const & p1)
    { return VectorNT((1.0 - coef) * p0 + coef * p1); }

}

#include "GeomXEdgeCalculations_t.hpp"

#endif
