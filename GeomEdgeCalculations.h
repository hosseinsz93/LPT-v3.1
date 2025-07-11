#ifndef _MKCOREEDGECALCULATIONS_H
#define _MKCOREEDGECALCULATIONS_H

#include <cmath> // for sqrt

#include "GeomVector.h"


namespace Geom
{
  class Id;
//   class TopologyData;
//   class GeometryData;
  class ConstEdgeVertices;

  // get the endpoint vertices of an edge
//   void edgeVertices(TopologyData *topology, Id const &edge,
//                                     Id &v1, Id &v2);

  // get the endpoint coordinates of an edge
//   void edgeCoordinates(GeometryData *geometry, TopologyData *topology, 
//                                       Id const &edge, Vector3d &coord1, 
//                                       Vector3d &coord2);

  // distance squared from edge to point
//   double edgePointDistanceSquared(GeometryData *geometry,
//                                   TopologyData *topology,
//                                   Id const &edge,
//                                   Vector3d const &querypt);
//   double edgePointDistanceSquared(GeometryData *geometry,
//                                   ConstEdgeVertices const &ev,
//                                   Vector3d const &querypt);
//   double edgePointDistanceSquared(GeometryData *geometry,
//                                   Id const &v1, Id const &v2,
//                                   Vector3d const &querypt);
  double edgePointDistanceSquared(Vector3d const &coord1,
                                  Vector3d const &coord2,
                                  Vector3d const &querypt,
                                  Vector3d *closestpoint = 0);

  // distance from edge to point
//   inline double edgePointDistance(GeometryData *geometry,
//                                   TopologyData *topology,
//                                   Id const &edge,
//                                   Vector3d const &querypt) {
//     return sqrt(edgePointDistanceSquared(geometry, topology, 
//                                          edge, querypt));
//   }
//   inline double edgePointDistance(GeometryData *geometry,
//                                   ConstEdgeVertices const &ev,
//                                   Vector3d const &querypt) {
//     return sqrt(edgePointDistanceSquared(geometry, ev, querypt));
//   }
//   inline double edgePointDistance(GeometryData *geometry,
//                                   Id const &v1, Id const &v2,
//                                   Vector3d const &querypt) {
//     return sqrt(edgePointDistanceSquared(geometry, v1, v2, querypt));
//   }
  inline double edgePointDistance(Vector3d const &coord1,
                                  Vector3d const &coord2,
                                  Vector3d const &querypt) {
    return sqrt(edgePointDistanceSquared(coord1, coord2, querypt));
  }

  // project point onto edge
  Vector3d edgeProjectPoint(Vector3d const &coord1,
                            Vector3d const &coord2,
                            Vector3d const &coord);

  // find ratio of projected point in terms of its distance from coord1.
  double edgeProjectRatio(Vector3d const &coord1,
                          Vector3d const &coord2,
                          Vector3d const &coord);

  double edgeProjectRatioExtended(Vector3d const &coord1,
                                  Vector3d const &coord2,
                                  Vector3d const &coord);

  Vector3d edgeCenter(Vector3d const &A, Vector3d const &B);

  Vector3d edgeDirectional(Vector3d const &A, Vector3d const &B);
}

#endif


// Automatic setting of emacs local variables.
// Local Variables:
// mode: C++
// tab-width: 8
// End:
