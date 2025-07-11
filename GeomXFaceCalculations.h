#ifndef _GEOMX_FACECALCULATIONS_H_
#define _GEOMX_FACECALCULATIONS_H_


#include <vector>

#include "GeomVector.h"

#include "GeomFaceCalculations.h"
#include "GeomId.h"
// #include "GeomXVertexCoordInterface.h"


namespace Geom
{
  class FaceVertices;
  class TopologyData;
}

namespace GeomX
{
  class VertexCoordInterface;

#if 0
  void faceCoordinates(VertexCoordInterface const * geometry,
                       GeomX::TopologyData const * topology,
                       GeomX::Id const &face,
                       std::vector< Vector3d > &coords);
#endif
  
  void getFaceCoords(VertexCoordInterface const * geom,
                     Geom::FaceVertices const & fv,
                     std::vector < Vector3d > & coords);

  Vector3d faceArea(VertexCoordInterface const * geom,
                            Geom::FaceVertices const & fv);

  Vector3d triangleAreaCoords(Vector3d const & p,
                              Vector3d const & c0,
                              Vector3d const & c1,
                              Vector3d const & c2);

  Vector3d triangleAreaCoords(Vector3d const & p,
                              Vector3d * cs);

  Vector3d triangleAreaCoords2D(Vector3d const & p,
                                Vector3d * cs);

  void triangleAreaCoordsConstrain(Vector3d const & p,
                                   Vector3d const & c0,
                                   Vector3d const & c1,
                                   Vector3d const & c2,
                                   Vector3d       & areaCoord);

  void triangleAreaCoordsConstrain(Vector3d const & p,
                                   Vector3d c[],
                                   Vector3d       & areaCoord);


  enum Sidedness
  {
    INSIDE   = 1<<0, // 0x01
    BOUNDARY = 1<<1, // 0x02
    OUTSIDE  = 1<<2, // 0x04
    UNKNOWN  = 1<<3, // 0x08
  };

  inline char const * toChar(Sidedness s)
  {
    switch (s)
    {
      case INSIDE:   return "INSIDE";
      case BOUNDARY: return "BOUNDARY";
      case OUTSIDE:  return "OUTSIDE";
      case UNKNOWN:  return "UNKNOWN";
    }
    return "Illegal Value!";
  }


  // Is the point inside or outside the given face (polygon).  This routine
  // works in 2D since this operation only makes sense in that system.
  // Dimensions higher than 2 are ignored.

  // C is a container class that contains VectorNT.  It has an operator[]
  // to access elements of the array and a size() procedure that returns
  // how many elements are in the array.  N must be 2 or greater.  T
  // must be float or double.
  template < typename C, int N, typename T >
  Sidedness pointSideofFace(C const & argPoly,
                            VectorNT const & test, double linearTol = 0.0);

  // Is the point inside or on the boundary of the polygon described by coords.

  // C is a container class that contains VectorNT.  It has an operator[]
  // to access elements of the array and a size() procedure that returns
  // how many elements are in the array.  N must be 3 and T must be double.
  template < typename C >
  bool pointInsideOnFace(C const & coords, Vector3d p);


  // Find the distance squared between a piont and a face.  The exact distance
  // is calculated.  This routine is for 2D only and dimensions higher than 2
  // are ignored.
  template < int N, typename T >
  double pointFaceDistance2(VectorNT p,
                            std::vector < VectorNT > const & fcoords,
                            double linearTol = 0.0);


  // This finds the almost actual distance between faces.  This is more exact
  // than the previous routine but it still makes a simplication.  The faces are
  // projected onto a plane defined by the face area vector and its centroid.

  // This is used by the interface routine.
  double faceFaceDistance(std::vector < Vector3d > const & coords0,
                          std::vector < Vector3d > const & coords1,
                          double tol, double * maxFacesDistance = 0);


  // See if the face is self-intersecting.
  bool faceIsSelfIntersecting(std::vector < Vector3d > const & cs);
}

#include "GeomXFaceCalculations_t.hpp"

#endif
