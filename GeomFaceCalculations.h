#ifndef _GEOM_FACECALCULATIONS_H
#define _GEOM_FACECALCULATIONS_H

#include <vector>
#include <cmath> // for sqrt
#include <limits>

#include "GeomUtility.h" // need it for ALMOST_ZERO definition
#include "GeomVector.h"



namespace Geom
{
  class Id;

  template <int N, typename T> class Vector;

  // Triangle access utility functions
  // void triangleVertices(TopologyData const * topology, Id const &face,
  //                       Id &v1, Id &v2, Id &v3);
  // void triangleCoordinates(GeometryData const * geometry,
  //                          TopologyData const * topology,
  //                          Id const &face, Vector3d &coord1,
  //                          Vector3d &coord2, Vector3d &coord3);

  // Triangle Area
  // To get magnitude or normalized area vector, use Vector3d.mag/.normalize

  // Vector3d triangleArea(GeometryData const * geometry,
  //                         TopologyData const * topology,
  //                         Id const &face);
  // Vector3d triangleArea(GeometryData const * geometry,
  //                         FaceVertices const &fv);
  // Vector3d triangleArea(GeometryData const *geometry,
  //                         Id const &v1, Id const &v2, Id const &v3);
  Vector3d triangleArea(Vector3d const &coord1,
                        Vector3d const &coord2,
                        Vector3d const &coord3);

  double triangleArea2D(Vector3d const &coord1,
                        Vector3d const &coord2,
                        Vector3d const &coord3);

  // Triangle Centroid

  // Vector3d triangleCentroid(GeometryData const * geometry,
  //                             TopologyData const * topology,
  //                             Id const &face);
  // Vector3d triangleCentroid(GeometryData const * geometry,
  //                             FaceVertices const &fv);
  // Vector3d triangleCentroid(GeometryData const * geometry,
  //                             Id const &v1, Id const &v2, Id const &v3);
  Vector3d triangleCentroid(Vector3d const &coord1,
                            Vector3d const &coord2,
                            Vector3d const &coord3);

  double trianglePerimeter(Vector3d const &A, 
                           Vector3d const &B, 
                           Vector3d const &C);

  Vector3d triangleNormal(Vector3d const &A, 
                          Vector3d const &B, 
                          Vector3d const &C); 

  Vector3d triangleIncenter(Vector3d const &A, 
                            Vector3d const &B, 
                            Vector3d const &C);

//   Vector3d triangleCircumcenter(Vector3d const &A,
//                                 Vector3d const &B, 
//                                 Vector3d const &C);

//   double triangleCircumradius(Vector3d const &A, 
//                               Vector3d const &B, 
//                               Vector3d const &C);

//  double triangleDistortionFactor(Vector3d const &A, 
//                                  Vector3d const &B, 
//                                  Vector3d const &C); 

  double triangleMeanRatio(Vector3d const &A, 
                           Vector3d const &B, 
                           Vector3d const &C);

 // Triangle Area and Centroid

  // void triangleAreaCentroid(GeometryData const * geometry,
  //                           TopologyData const * topology,
  //                           Id const &face,
  //                           Vector3d &area, Vector3d &centroid);
  // void triangleAreaCentroid(GeometryData *geometry, FaceVertices const &fv,
  //                           Vector3d &area, Vector3d &centroid);
  // void triangleAreaCentroid(GeometryData const * geometry,
  //                           Id const &v1, Id const &v2, Id const &v3,
  //                           Vector3d &area, Vector3d &centroid);
  void triangleAreaCentroid(Vector3d const &coord1,
                            Vector3d const &coord2,
                            Vector3d const &coord3,
                            Vector3d &area, Vector3d &centroid);


  // Triangle quality calculations

  double triangleQuality1(Vector3d const &a,
                          Vector3d const &b, Vector3d const &c);



  // squared closest distance from triangle to point 
  // double trianglePointDistanceSquared(GeometryData *geometry,
  //                                     FaceVertices const &fv,
  //                                     Vector3d const &querypt);

  double trianglePointDistanceSquared(Vector3d const &coord1,
				      Vector3d const &coord2,
				      Vector3d const &coord3,
				      Vector3d const &querypt);

  double trianglePointDistanceSquared(Vector3d const &coord1,
                                      Vector3d const &coord2,
                                      Vector3d const &coord3,
                                      Vector3d const &querypt,
                                      Vector3d *ptOnTriangle);
  // closest distance from triangle to point
  // inline double trianglePointDistance(GeometryData *geometry,
  //                                     FaceVertices const &fv,
  //                                     Vector3d const &querypt) {
  //   return sqrt(trianglePointDistanceSquared(geometry, fv, querypt));
  // }
  inline double trianglePointDistance(Vector3d const &coord1,
                                      Vector3d const &coord2,
                                      Vector3d const &coord3,
                                      Vector3d const &querypt) {
    return sqrt(trianglePointDistanceSquared(coord1, coord2,
                                             coord3, querypt));
  }
    
  // baricentric coordinates
  Vector3d computeBaricentricCoordinatesFromReal(Vector3d const &coord0,
                                                 Vector3d const &coord1,
                                                 Vector3d const &coord2,
                                                 Vector3d const &coord);
  Vector3d computeRealCoordinatesFromBaricentric(Vector3d const &coord0,
                                                 Vector3d const &coord1,
                                                 Vector3d const &coord2,
                                                 Vector3d const &baric);
  // trilinear coordinates
  Vector3d computeTrilinearCoordinatesFromReal(Vector3d const &coord0,
                                               Vector3d const &coord1,
                                               Vector3d const &coord2,
                                               Vector3d const &coord);
  Vector3d computeRealCoordinatesFromTrilinear(Vector3d const &coord0,
                                               Vector3d const &coord1,
                                               Vector3d const &coord2,
                                               Vector3d const &trilin);
  Vector3d computeRealCoordinatesFromBlendedBaricentricTrilinear(
                                               Vector3d const &coord0,
					       Vector3d const &coord1,
					       Vector3d const &coord2,
					       Vector3d const &baric,
                                               Vector3d const &trilin);

  // Face access utility functions

  // void faceVertices(TopologyData *topology, Id const &face,
  //                   std::vector<Id> &vertices);
  // void faceCoordinates(GeometryData const * geometry,
  //                      TopologyData const * topology,
  //                      Id const &face, std::vector< Vector3d > &coord);




  // Face Area
  // To get magnitude or normalized area vector, use Vector3d.mag/.normalize

  // Vector3d faceArea(GeometryData *geometry,
  //                                    TopologyData *topology, Id const &face);
  // Vector3d faceArea(GeometryData *geometry, std::vector<Id> &vertices);
  // Vector3d faceArea(GeometryData const * geometry,
  //                                          FaceVertices const &fv);
  Vector3d faceArea(std::vector<Vector3d > const &coord);


  // Face Centroid

  // Vector3d faceCentroid(GeometryData *geometry, TopologyData *topology,
  //                         Id const &face, double AreaTol=ALMOST_ZERO);
  // Vector3d faceCentroid(GeometryData *geometry, std::vector<Id> &vertices,
  //                         double AreaTol=ALMOST_ZERO);
  // Vector3d faceCentroid(GeometryData *geometry, FaceVertices const &fv,
  //                         double AreaTol=ALMOST_ZERO);
  Vector3d faceCentroid(std::vector<Vector3d > const &coord,
                          double AreaTol=ALMOST_ZERO);




  // Face Area and Centroid

  // void faceAreaCentroid(GeometryData *geometry, TopologyData *topology,
  //                       Id const &face,
  //                       Vector3d &area, Vector3d &centroid,
  //                       double AreaTol=ALMOST_ZERO);
  // void faceAreaCentroid(GeometryData *geometry, std::vector<Id> &vertices,
  //                       Vector3d &area, Vector3d &centroid,
  //                       double AreaTol=ALMOST_ZERO);
  // void faceAreaCentroid(GeometryData *geometry, FaceVertices const &fv,
  //                       Vector3d &area, Vector3d &centroid,
  //                       double AreaTol=ALMOST_ZERO);
  void faceAreaCentroid(std::vector<Vector3d > const &coord,
                        Vector3d &area, Vector3d &centroid,
                        double AreaTol=ALMOST_ZERO);


  // squared closest distance from face to point
  // double facePointDistanceSquared(GeometryData *geometry, 
  //                                                TopologyData *topology,
  //                                                Id const &face,
  //                                                Vector3d const &querypt);
  // double facePointDistanceSquared(GeometryData *geometry, 
  //                                                FaceVertices const &fv,
  //                                                Vector3d const &querypt);
  double facePointDistanceSquared(std::vector<Vector3d > const &coord,
                                  Vector3d const &querypt);

  // closest distance from face to point
  // inline double facePointDistance(GeometryData *geometry, 
  //                                 TopologyData *topology,
  //                                 Id const &face,
  //                                 Vector3d const &querypt) {
  //   return sqrt(facePointDistanceSquared(geometry, topology, face, 
  //                                        querypt));
  // }
    
  // inline double facePointDistance(GeometryData *geometry, 
  //                                 FaceVertices const &fv,
  //                                 Vector3d const &querypt) {
  //   return sqrt(facePointDistanceSquared(geometry, fv, querypt));
  // }
    
  inline double facePointDistance(std::vector<Vector3d > const &coord,
                                  Vector3d const &querypt) {
    return sqrt(facePointDistanceSquared(coord, querypt));
  }


  // face planarity
  // double facePlanarity(GeometryData *geometry,
  // 		      TopologyData *topology, Id const &face);
  // double facePlanarity(GeometryData *geometry, std::vector<Id> &vertices);
  // double facePlanarity(GeometryData const * geometry,
  // 		      FaceVertices const &fv);
  double facePlanarity(std::vector<Vector3d > const &coord);
  double facePlanarity(std::vector<Vector3d > const &coord,
                       Vector3d const &faceArea,
                       Vector3d const &faceCentroid);

  // face twist
  // double faceTwist(GeometryData *geometry,
  // 		  TopologyData *topology, Id const &face);
  // double faceTwist(GeometryData *geometry, std::vector<Id> &vertices);
  // double faceTwist(GeometryData const * geometry,
  // 		  FaceVertices const &fv);
  double faceTwist(std::vector<Vector3d > const &coord);
  double faceTwist(std::vector<Vector3d > const &coord,
                   Vector3d const &faceArea,
                   Vector3d const &faceCentroid);


  // Calculate the minimum edge length for a face.
  // Vec requires a subtract and length (magnitude) function.
  template < int N, class T, template < int O, class U > class Vec >
  T faceMinEdgeLeng(Vec<N,T> const * vec, int pn)
  {
    if (pn <= 1)
      return 0.0;

    int i, p;
    T minLeng = std::numeric_limits < T >::max();
    for (i=0, p=pn-1; i<pn; p=i, ++i)
    {
      T leng = Vec<N,T>(vec[p] - vec[i]).length();
      minLeng = std::min(leng, minLeng);
    }
    return minLeng;
  }
}

#endif


// Automatic setting of emacs local variables.
// Local Variables:
// mode: C++
// tab-width: 8
// End:
