// This file and its contents are the property of adapco, Ltd. and contain
// confidential and proprietary information.  The unauthorized use,
// distribution, or duplication of this information is prohibited.
// All rights reserved.
//
// Copyright (C) 2004 by adapco, Ltd.

#include "GeomFaceCalculations.h"

#include <limits>

// Temporary change to make parallel trimming abort, not hang. -Larry
// #include "neo/Assert.h"

#include "GeomEdgeCalculations.h"

// #include "common/Tensor.h"

#include "GeomException.h"
#include "GeomUtility.h"

#include "GeomTriangle3.h"
#include "GeomPoint3.h"

#include "GeomId.h"
// #include "GeometryData.h"
// #include "TopologyDefs.h"
// #include "TopologyData.h"

namespace Geom
{

  // Triangle access utility functions

  /* Get the vertices of a triangle.  This routine loops through the different
   * topology pointers to extract the vertices of the triangular face 'face'.
   * The routine throws an exception if:
   * 1) Face has less than 3 vertices/edges
   * 2) Face has greater than 3 vertices/edges
   * 3) Insufficient/incorrect information
   */

  // void triangleVertices(TopologyData const * topology, Id const &face,
  //                       Id &v1, Id &v2, Id &v3)
  // {
  //   if (topology->existsFtoV(face))
  //     {
  //       ConstFaceVertices fv = topology->getFtoV(face);
  //       if (fv.size() < 3)
  //         throw Exception("triangleVertices:  less than 3 vertices");
  //       else if (fv.size() == 3)
  //         {
  //           v1 = fv[0];
  //           v2 = fv[1];
  //           v3 = fv[2];
  //           return;
  //         }
  //       else
  //         throw Exception("triangleVertices:  more than 3 vertices");
  //     }
  //   else if (topology->existsFtoE(face))
  //     {
  //       ConstFaceEdges fe = topology->getFtoE(face);
  //       if (fe.size() < 3)
  //         throw Exception("triangleVertices:  less than 3 edges");
  //       else if (fe.size() == 3)
  //         {
  //           if (!topology->existsEtoV(fe[0]) ||
  //              !topology->existsEtoV(fe[1]) ||
  //              !topology->existsEtoV(fe[2]))
  //             throw Exception("triangleVertices:  insufficient information");
  //           ConstEdgeVertices e1 = topology->getEtoV(fe[0]);
  //           ConstEdgeVertices e2 = topology->getEtoV(fe[1]);
  //           // Get v1 (the common vertex between e1 and e2)
  //           if (e1[0] == e2[0] || e1[0] == e2[1])
  //             v1 = e1[0];
  //           else if (e1[1] == e2[0] || e1[1] == e2[1])
  //             v1 = e1[1];
  //           else
  //             throw Exception("triangleVertices:  incorrect information");
  //           // Get v2 and v3 (assumption: face is truely closed!)
  //           v2 = Id(e2[0].toUnsigned() + e2[1].toUnsigned() - v1.toUnsigned()); // other endpoint of e2
  //           v3 = Id(e1[0].toUnsigned() + e1[1].toUnsigned() - v1.toUnsigned()); // other endpoint of e1 (no need for e3)
  //         }
  //       else
  //         throw Exception("triangleVertices:  more than 3 edges");
  //     }
  //   else
  //     throw Exception("triangleVertices:  insufficient information");

  //    return;
  // }




  // void triangleCoordinates(GeometryData const * geometry,
  //                          TopologyData const * topology,
  //                          Id const &face, Vector3d &coord1,
  //                          Vector3d &coord2, Vector3d &coord3)
  // {
  //   Id v1, v2, v3;
  //   triangleVertices(topology, face, v1, v2, v3);
  //   coord1 = geometry->getPoint(v1);
  //   coord2 = geometry->getPoint(v2);
  //   coord3 = geometry->getPoint(v3);
  //   return;
  // }




  // Triangle Area

  // Vector3d triangleArea(GeometryData const * geometry,
  //                         TopologyData const * topology,
  //                         Id const &face)
  // {
  //   Id v1, v2, v3;

  //    triangleVertices(topology, face, v1, v2, v3);
  //   return triangleArea(geometry, v1, v2, v3);
  // }




  // Vector3d triangleArea(GeometryData const * geometry,
  //                         FaceVertices const &fv)
  // {
  //   if (fv.size() != 3)
  //     throw Exception("triangleArea:  not a triangular face");

  //    return triangleArea(geometry, fv[0], fv[1], fv[2]);
  // }




  // Vector3d triangleArea(GeometryData const * geometry,
  //                         Id const &v1, Id const &v2, Id const &v3)
  // {
  //   return triangleArea(geometry->getPoint(v1),
  //                       geometry->getPoint(v2),
  //                       geometry->getPoint(v3));
  // }




  Vector3d triangleArea(Vector3d const &coord1,
                        Vector3d const &coord2,
                        Vector3d const &coord3)
  {
    Vector3d vec1(coord3 - coord2);
    Vector3d vec2(coord1 - coord2);

    Vector3d area = vec1.cross(vec2);
    return area*=0.5;
  }

  
  double triangleArea2D(Vector3d const &coord1,
                        Vector3d const &coord2,
                        Vector3d const &coord3)
  {
    Vector3d vec1(coord3 - coord2); vec1(2) = 0.0;
    Vector3d vec2(coord1 - coord2); vec2(2) = 0.0;

    Vector3d area = vec1.cross(vec2);
    return 0.5 * area(2);
  }

  // Triangle Centroid

  // Vector3d triangleCentroid(GeometryData const * geometry,
  //                             TopologyData const * topology,
  //                             Id const &face)
  // {
  //   Id v1, v2, v3;

  //    triangleVertices(topology, face, v1, v2, v3);
  //   return triangleCentroid(geometry, v1, v2, v3);
  // }




  // Vector3d triangleCentroid(GeometryData const * geometry,
  //                             FaceVertices const &fv)
  // {
  //   if (fv.size() != 3)
  //     throw Exception("triangleCentroid:  not a triangular face");

  //    return triangleCentroid(geometry, fv[0], fv[1], fv[2]);
  // }




  // Vector3d triangleCentroid(GeometryData const * geometry,
  //                             Id const &v1, Id const &v2, Id const &v3)
  // {
  //   return triangleCentroid(geometry->getPoint(v1),
  //                           geometry->getPoint(v2),
  //                           geometry->getPoint(v3));
  // }




  Vector3d triangleCentroid(Vector3d const &coord1,
                            Vector3d const &coord2,
                            Vector3d const &coord3)
  {
    Vector3d centroid(coord1 + coord2 + coord3);
    return centroid *= (1.0/3.0);
  }

  // Triangle area and Centroid

  // void triangleAreaCentroid(GeometryData const * geometry,
  //                           TopologyData const * topology,
  //                           Id const &face,
  //                           Vector3d &area, Vector3d &centroid)
  // {
  //   Id v1, v2, v3;

  //    triangleVertices(topology, face, v1, v2, v3);
  //   triangleAreaCentroid(geometry, v1, v2, v3, area, centroid);
  //   return;
  // }

  // void triangleAreaCentroid(GeometryData *geometry, FaceVertices const &fv,
  //                           Vector3d &area, Vector3d &centroid)
  // {
  //   if (fv.size() != 3)
  //     throw Exception("triangleAreaCentroid:  not a triangular face");

  //    triangleAreaCentroid(geometry, fv[0], fv[1], fv[2], area, centroid);
  //   return;
  // }

  // void triangleAreaCentroid(GeometryData const * geometry,
  //                           Id const &v1, Id const &v2, Id const &v3,
  //                           Vector3d &area, Vector3d &centroid)
  // {
  //   triangleAreaCentroid(geometry->getPoint(v1),
  //                        geometry->getPoint(v2),
  //                        geometry->getPoint(v3),
  //                        area, centroid);
  //   return;
  // }

  void triangleAreaCentroid(Vector3d const &coord1,
                            Vector3d const &coord2,
                            Vector3d const &coord3,
                            Vector3d &area, Vector3d &centroid)
  {
    area     = triangleArea(coord1, coord2, coord3);
    centroid = triangleCentroid(coord1, coord2, coord3);
    return;
  }

  double trianglePerimeter(Vector3d const &A, 
                           Vector3d const &B, 
                           Vector3d const &C) 
  {  
    return A.distance(B) 
         + B.distance(C) 
         + C.distance(A);
  }

  Vector3d triangleNormal(Vector3d const &A, 
                          Vector3d const &B, 
                          Vector3d const &C) 
  {  
    Vector3d u(B - A);
    Vector3d v(C - A);
    Vector3d uv = u.cross(v);
    uv.normalize();
    return uv;
  }

  Vector3d triangleIncenter(Vector3d const &A, 
                            Vector3d const &B, 
                            Vector3d const &C) 
  {
    double d1 = B.distance(C);
    double d2 = C.distance(A);
    double d3 = A.distance(B);
    return Vector3d((A * d1 + B * d2 + C * d3) / (d1 + d2 + d3));
  }
  
//   Vector3d triangleCircumcenter(Vector3d const &A,
//                                         Vector3d const &B, 
//                                         Vector3d const &C) 
//   {
//     Tensor<3,double> a;
//     Vector3d b, rc;
// 
//     a(0,0) = A(0) - C(0);
//     a(0,1) = A(1) - C(1);
//     a(0,2) = A(2) - C(2);
//     a(1,0) = B(0) - C(0);
//     a(1,1) = B(1) - C(1);
//     a(1,2) = B(2) - C(2);
//     a(2,0) = a(0,1) * a(1,2) - a(0,2) * a(1,1);
//     a(2,1) = a(0,2) * a(1,0) - a(0,0) * a(1,2);
//     a(2,2) = a(0,0) * a(1,1) - a(0,1) * a(1,0);
// 
//     b(0) = (a(0,0) * (A(0) + C(0)) + a(0,1) * (A(1) + C(1)) + a(0,2) * (A(2) + C(2))) / 2;
//     b(1) = (a(1,0) * (B(0) + C(0)) + a(1,1) * (B(1) + C(1)) + a(1,2) * (B(2) + C(2))) / 2;
//     b(2) =  a(2,0) *         C(0)  + a(2,1)         * C(1)  + a(2,2)         * C(2);
// 
//     double det = a.determinant();
//     //if (fsbs(det) < 1.0e-30) return triangleCentroid(A, B, C);    
//     if (det == 0.0) return triangleCentroid(A, B, C);    
//      
//     rc = Tensor<3,double>(a.adjoint()/det).dot(b);
//     return rc;
//   }

//   double triangleCircumradius(Vector3d const &A, 
//                               Vector3d const &B, 
//                               Vector3d const &C) 
//   {
//     return triangleCircumcenter(A, B, C).distance(A);
//   }

//  double triangleDistortionFactor(Vector3d const &A, 
//                                  Vector3d const &B, 
//                                  Vector3d const &C) 
//  {
//    double area = triangleArea(A, B, C).length();
//    if (area == 0.0) return 0.0;
//    double radius = triangleCircumradius(A, B, C);
//    if (radius == 0.0) return 0.0;
//    double const factor = 0.769800358919501;  // 4/3/3^(1/2)
//    double q = factor * area / (radius * radius);
//    if (q > 1.0) q = 1.0;
//    return q;
//  }  

  double triangleMeanRatio(Vector3d const &A, 
                           Vector3d const &B, 
                           Vector3d const &C)
  {
    double area = triangleArea(A, B, C).length();
    if (area == 0.0) return 0.0;
    double dd = A.distanceSquared(B)
              + B.distanceSquared(C)
              + C.distanceSquared(A);
    if (dd == 0.0) return 0.0;
    double const factor = 6.928203230275509;  // 12/3^(1/2)
    double q = factor * area / dd;
    if (q > 1.0) q = 1.0;
    return q;
  }  

  // Triangle quality calculations

  double triangleQuality1(Vector3d const &a,
                          Vector3d const &b, Vector3d const &c)
  {
    return triangleMeanRatio(a, b, c);
  }


  // squared closest distance from triangle to point 

  // double trianglePointDistanceSquared(GeometryData *geometry,
  //                                     FaceVertices const &fv,
  //                                     Vector3d const &querypt)
  // {
  //   if (fv.size() != 3)
  //     throw Exception(
  //                     "trianglePointDistanceSquared:  not a triangular face");
  //   return trianglePointDistanceSquared(geometry->at(fv[0]),
  //                                       geometry->at(fv[1]),
  //                                       geometry->at(fv[2]),
  //                                       querypt);
  // }

  double trianglePointDistanceSquared(Vector3d const &coord1,
                                      Vector3d const &coord2,
                                      Vector3d const &coord3,
                                      Vector3d const &querypt)
  {   
    return trianglePointDistanceSquared(coord1,coord2,coord3,querypt,0);
  }

  double trianglePointDistanceSquared(Vector3d const &coord1,
                                      Vector3d const &coord2,
                                      Vector3d const &coord3,
                                      Vector3d const &querypt,
                                      Vector3d *ptOnTriangle)
  {

    // algorithm lifted from Wild Magick Library 
    // (www.magick-software.com)
    Vector3d edge0(coord2 - coord1);
    Vector3d edge1(coord3 - coord1);

    Vector3d kDiff(coord1 - querypt);
    double fA00 = edge0.mag2();
    double fA01 = edge0.dot(edge1);
    double fA11 = edge1.mag2();
    double fB0 = kDiff.dot(edge0);
    double fB1 = kDiff.dot(edge1);
    double fC = kDiff.mag2();
    double fDet = fabs(fA00*fA11-fA01*fA01);

    if (fDet < ALMOST_ZERO) {
      // triangle is degenerate -- either all three points are coincient
      // or colinear.  

      // find longest edge
      std::vector<Vector3d const *> v(3);
      v[0] = &coord1;
      v[1] = &coord2;
      v[2] = &coord3;

      double maxEdgeLengthSq = 0.0;
      int longestEdge = 0;
      for (int i = 0; i<3; ++i) {
	double edgeLengthSq = (*v[i]).distanceSquared(*v[(i+1)%3]); 
	if (edgeLengthSq > maxEdgeLengthSq) {
	  maxEdgeLengthSq = edgeLengthSq;
	  longestEdge = i;
	}
      }

      return edgePointDistanceSquared(*v[longestEdge], *v[(longestEdge+1)%3],
				      querypt,ptOnTriangle);

    }


    double fS = fA01*fB1-fA11*fB0;
    double fT = fA01*fB0-fA00*fB1;
    double fSqrDist;

    if ( fS + fT <= fDet )
      {
        if ( fS < 0.0 )
          {
            if ( fT < 0.0 )  // region 4
              {
                if ( fB0 < 0.0 )
                  {
                    fT = 0.0;
                    if ( -fB0 >= fA00 )
                      {
                        fS = 1.0;
                        fSqrDist = fA00+2.0*fB0+fC;
                      }
                    else
                      {
                        fS = -fB0/fA00;
                        fSqrDist = fB0*fS+fC;
                      }
                  }
                else
                  {
                    fS = 0.0;
                    if ( fB1 >= 0.0 )
                      {
                        fT = 0.0;
                        fSqrDist = fC;
                      }
                    else if ( -fB1 >= fA11 )
                      {
                        fT = 1.0;
                        fSqrDist = fA11+2.0*fB1+fC;
                      }
                    else
                      {
                        fT = -fB1/fA11;
                        fSqrDist = fB1*fT+fC;
                      }
                  }
              }
            else  // region 3
              {
                fS = 0.0;
                if ( fB1 >= 0.0 )
                  {
                    fT = 0.0;
                    fSqrDist = fC;
                  }
                else if ( -fB1 >= fA11 )
                  {
                    fT = 1.0;
                    fSqrDist = fA11+2.0*fB1+fC;
                  }
                else
                  {
                    fT = -fB1/fA11;
                    fSqrDist = fB1*fT+fC;
                  }
              }
          }
        else if ( fT < 0.0 )  // region 5
          {
            fT = 0.0;
            if ( fB0 >= 0.0 )
              {
                fS = 0.0;
                fSqrDist = fC;
              }
            else if ( -fB0 >= fA00 )
              {
                fS = 1.0;
                fSqrDist = fA00+2.0*fB0+fC;
              }
            else
              {
                fS = -fB0/fA00;
                fSqrDist = fB0*fS+fC;
              }
          }
        else  // region 0
          {
            // minimum at interior point
            double fInvDet = 1.0/fDet;
            fS *= fInvDet;
            fT *= fInvDet;
            fSqrDist = fS*(fA00*fS+fA01*fT+2.0*fB0) +
              fT*(fA01*fS+fA11*fT+2.0*fB1)+fC;
          }
      }
    else
      {
        double fTmp0, fTmp1, fNumer, fDenom;
            
        if ( fS < 0.0 )  // region 2
          {
            fTmp0 = fA01 + fB0;
            fTmp1 = fA11 + fB1;
            if ( fTmp1 > fTmp0 )
              {
                fNumer = fTmp1 - fTmp0;
                fDenom = fA00-2.0f*fA01+fA11;
                if ( fNumer >= fDenom )
                  {
                    fS = 1.0;
                    fT = 0.0;
                    fSqrDist = fA00+2.0*fB0+fC;
                  }
                else
                  {
                    fS = fNumer/fDenom;
                    fT = 1.0 - fS;
                    fSqrDist = fS*(fA00*fS+fA01*fT+2.0f*fB0) +
                      fT*(fA01*fS+fA11*fT+2.0*fB1)+fC;
                  }
              }
            else
              {
                fS = 0.0;
                if ( fTmp1 <= 0.0 )
                  {
                    fT = 1.0;
                    fSqrDist = fA11+2.0*fB1+fC;
                  }
                else if ( fB1 >= 0.0 )
                  {
                    fT = 0.0;
                    fSqrDist = fC;
                  }
                else
                  {
                    fT = -fB1/fA11;
                    fSqrDist = fB1*fT+fC;
                  }
              }
          }
        else if ( fT < 0.0 )  // region 6
          {
            fTmp0 = fA01 + fB1;
            fTmp1 = fA00 + fB0;
            if ( fTmp1 > fTmp0 )
              {
                fNumer = fTmp1 - fTmp0;
                fDenom = fA00-2.0*fA01+fA11;
                if ( fNumer >= fDenom )
                  {
                    fT = 1.0;
                    fS = 0.0;
                    fSqrDist = fA11+2.0*fB1+fC;
                  }
                else
                  {
                    fT = fNumer/fDenom;
                    fS = 1.0 - fT;
                    fSqrDist = fS*(fA00*fS+fA01*fT+2.0*fB0) +
                      fT*(fA01*fS+fA11*fT+2.0*fB1)+fC;
                  }
              }
            else
              {
                fT = 0.0;
                if ( fTmp1 <= 0.0 )
                  {
                    fS = 1.0;
                    fSqrDist = fA00+2.0*fB0+fC;
                  }
                else if ( fB0 >= 0.0 )
                  {
                    fS = 0.0;
                    fSqrDist = fC;
                  }
                else
                  {
                    fS = -fB0/fA00;
                    fSqrDist = fB0*fS+fC;
                  }
              }
          }
        else  // region 1
          {
            fNumer = fA11 + fB1 - fA01 - fB0;
            if ( fNumer <= 0.0 )
              {
                fS = 0.0;
                fT = 1.0;
                fSqrDist = fA11+2.0*fB1+fC;
              }
            else
              {
                fDenom = fA00-2.0f*fA01+fA11;
                if ( fNumer >= fDenom )
                  {
                    fS = 1.0;
                    fT = 0.0;
                    fSqrDist = fA00+2.0*fB0+fC;
                  }
                else
                  {
                    fS = fNumer/fDenom;
                    fT = 1.0 - fS;
                    fSqrDist = fS*(fA00*fS+fA01*fT+2.0*fB0) +
                      fT*(fA01*fS+fA11*fT+2.0*fB1)+fC;
                  }
              }
          }
      }
        
    if (ptOnTriangle)
      {
        *ptOnTriangle = coord1 + edge0*fS + edge1*fT;
      }
    return fabs(fSqrDist);

  }



  // triangle coordinate systems (baricentric , trilinear )


  /* NOTE:  These routines allows the conversion to and from real, baricentric
            and trilinear

	    the trilinear coordinates are not 'exact' i.e. not normalized while
            the baricentric coordinates are normalized.
	    Each triplet of trilinear coordinates with the same relative ratios 
	    represents the same point in R3 (on the same triangle).

	    coord0,coord1,coord2 are the coordinates in real space of the 
            vertices of the triangle.
	    coord  are the coordinates of the point in real space
	    baric  are the baricentric coordinates
	    trilin are the trilinear coordinates
   */


  // baricentric coordinates
  Vector3d computeBaricentricCoordinatesFromReal(Vector3d const &coord0,
                                                 Vector3d const &coord1,
                                                 Vector3d const &coord2,
                                                 Vector3d const &coord)
  {
    Vector3d e0, e1, e2;
    Vector3d a0, a1, a2;
    double A0,A1,A2,AT,AINV;
    Vector3d b;

    // the interior edges
    e0=coord0-coord;
    e1=coord1-coord;
    e2=coord2-coord;

    // compute area vectors
    a0=e1.cross(e2);
    a1=e2.cross(e0);
    a2=e0.cross(e1);
    
    // evaluate areas
    A0=a0.length();
    A1=a1.length();
    A2=a2.length();
    AT=A0+A1+A2;

    // check if zero area
    if (AT==0.0) {
      Exception("Cannot compute baricentric coordinates on zero area triangle");
    }

    // compute baricentric coords
    AINV=1.0/AT;
    b(0)=A0*AINV;
    b(1)=A1*AINV;
    b(2)=A2*AINV;

    return b;
  }

  Vector3d computeRealCoordinatesFromBaricentric(Vector3d const &coord0,
                                                 Vector3d const &coord1,
                                                 Vector3d const &coord2,
                                                 Vector3d const &baric)
  {
    Vector3d r(coord0*baric(0) + coord1*baric(1) + coord2*baric(2));
    return r;
  }

  // trilinear coordinates
  Vector3d computeTrilinearCoordinatesFromReal(Vector3d const &coord0,
                                               Vector3d const &coord1,
                                               Vector3d const &coord2,
                                               Vector3d const &coord)
  {
    Vector3d e0, e1, e2;
    double l0,l1,l2;
    Vector3d trilin, baric;
    
    // the edge vectors
    e0=coord2-coord1;
    e1=coord0-coord2;
    e2=coord1-coord0;

    // the edge lengths
    l0=e0.length();
    l1=e1.length();
    l2=e2.length();

    // compute baricentric coords
    baric= computeBaricentricCoordinatesFromReal(coord0,coord1,coord2, coord);
    

    // compute trilinear coordinates
    trilin(0)=baric(0)/l0;
    trilin(1)=baric(1)/l1;
    trilin(2)=baric(2)/l2;
    
    return trilin;

  }

  Vector3d computeRealCoordinatesFromTrilinear(Vector3d const &coord0,
                                               Vector3d const &coord1,
                                               Vector3d const &coord2,
                                               Vector3d const &trilin)
  {
    Vector3d r, baric;
    double bt,binv;
    Vector3d e0, e1, e2;
    double l0,l1,l2;

    // the edge vectors
    e0=coord2-coord1;
    e1=coord0-coord2;
    e2=coord1-coord0;

    //the edge lengths
    l0=e0.length();
    l1=e1.length();
    l2=e2.length();

    //compute non normalized baricentric coordinates
    baric(0)=trilin(0)*l0;
    baric(1)=trilin(1)*l1;
    baric(2)=trilin(2)*l2;
    
    //normalize
    bt=baric(0)+baric(1)+baric(2);
    binv=1.0/bt;
    baric*=binv;

    //compute real from baricentric
    return computeRealCoordinatesFromBaricentric(coord0,coord1,coord2,baric);

  }

  Vector3d computeRealCoordinatesFromBlendedBaricentricTrilinear(
                                                        Vector3d const &coord0,
                                                        Vector3d const &coord1,
                                                        Vector3d const &coord2,
                                                        Vector3d const &baric,
                                                        Vector3d const &trilin)
  {
    // It blends the baricentric / trilinear coordinates in such a way that 
    // if the trilinear refers to the incenter (1,1,1 or any equal value triplet) 
    // then the coordinates of the incenter are returned.
    // If a point is on the boundary it returns the coordinates corresponding to the baricentric coordinates
    // It blends linearly inbetween.

    Vector3d r, baric_trl, baric_bln;
    double bt,binv;
    Vector3d e0, e1, e2;
    double l0,l1,l2;

    //compute blend factor
    double st=trilin(0)+trilin(1)+trilin(2);
    Vector3d norm_trilin(trilin*3.0/st);
    double blend_factor = norm_trilin(0)*norm_trilin(1)*norm_trilin(2);  // it is 1 at incenter and 0 at boundary

    //the edge vectors
    e0=coord2-coord1;
    e1=coord0-coord2;
    e2=coord1-coord0;

    //the edge lengths
    l0=e0.length();
    l1=e1.length();
    l2=e2.length();

    //compute non normalized baricentric coordinates corresponding to pure trilinear
    baric_trl(0)=trilin(0)*l0;
    baric_trl(1)=trilin(1)*l1;
    baric_trl(2)=trilin(2)*l2;
    
    //normalize
    bt=baric_trl(0)+baric_trl(1)+baric_trl(2);
    binv=1.0/bt;
    baric_trl*=binv;

    //blend
    baric_bln = baric_trl*blend_factor + baric*(1.0-blend_factor);

    //normalize
    bt=baric_bln(0)+baric_bln(1)+baric_bln(2);
    binv=1.0/bt;
    baric_bln*=binv;


    //compute real from baricentric
    return computeRealCoordinatesFromBaricentric(coord0,coord1,coord2,baric_bln);

  }




  // Face access utility functions

  /* Get the vertices of a face.  This routine loops through the different
   * topology pointers to extract the vertices of the face 'face'.
   *
   * NOTE:  The routines faceVertices and faceCoordinates are EXACTLY the same
   *        except for the fact that one returns vertex ids and one returns
   *        vertex coordinates.  The reason I have the code duplicated is that
   *        if all I need is vertex coordinates, there is no need to construct
   *        an stl vector composed of ids because I then need to construct
   *        another stl vector composed of coordinates.  If at some point we
   *        find that the construction of the stl vector is cheap, the
   *        faceCoordinates routine can simply call faceVertices and go from
   *        there...
   */

  // SEE NOTE ABOVE
  // void faceVertices(TopologyData *topology, Id const &face,
  //                   std::vector<Id> &vertices)
  // {
  //   vertices.clear(); // QUESTION: Is this expensive or simply setting size=0?

  //    if (topology->existsFtoV(face))
  //     {
  //       ConstFaceVertices fv = topology->getFtoV(face);
  //       vertices.reserve(fv.size());
  //       for (unsigned i=0; i<fv.size(); i++)
  //         vertices.push_back(fv[i]);
  //     }
  //   else if (topology->existsFtoE(face))
  //     {
  //       ConstFaceEdges fe = topology->getFtoE(face);
  //       for (unsigned i=0; i<fe.size(); i++)
  //         if (!topology->existsEtoV(fe[i]))
  //           throw Exception("faceVertices:  insufficient information");
  //       if (fe.size() == 0) // empty -- nothing to do
  //         ;
  //       else if (fe.size() == 1) // only one edge -- just get its vertices
  //         {
  //           ConstEdgeVertices e1 = topology->getEtoV(fe[0]);
  //           vertices.reserve(2);
  //           vertices.push_back(e1[0]);
  //           vertices.push_back(e1[1]);
  //         }
  //       else // Regular part
  //         {
  //           Id v1, v2, vn;
  //           ConstEdgeVertices e1 = topology->getEtoV(fe[0]);
  //           ConstEdgeVertices e2 = topology->getEtoV(fe[1]);
  //           // Get v1 (the common vertex between e1 and e2)
  //           if (e1[0] == e2[0] || e1[0] == e2[1])
  //             v1 = e1[0];
  //           else if (e1[1] == e2[0] || e1[1] == e2[1])
  //             v1 = e1[1];
  //           else
  //             throw Exception("faceVertices:  incorrect information");
  //           // Get the other two vertices of these edges
  //           v2 = Id(e2[0].toUnsigned() + e2[1].toUnsigned() - v1.toUnsigned()); // second vertex of face
  //           vn = Id(e1[0].toUnsigned() + e1[1].toUnsigned() - v1.toUnsigned()); // last vertex of face
  //           // Fill the vertices vector by traversing through face edges
  //           vertices.reserve(fe.size());
  //           vertices.push_back(v1);
  //           vertices.push_back(v2);
  //           for (unsigned i=2; i<fe.size()-1; i++)
  //             {
  //               ConstEdgeVertices e = topology->getEtoV(fe[i]);
  //               v2 = Id(e[0].toUnsigned() + e[1].toUnsigned() - v2.toUnsigned());
  //               vertices.push_back(v2);
  //             }
  //           vertices.push_back(vn);
  //         }
  //     }
  //   else
  //     throw Exception("faceVertices:  insufficient information");

  //    return;
  // }




  // SEE NOTE ABOVE
  // void faceCoordinates(GeometryData const * geometry,
  //                      TopologyData const * topology,
  //                      Id const &face, std::vector < Vector3d > &coord)
  // {
  //   coord.clear(); // QUESTION: Is this expensive or simply setting size=0?

  //    if (topology->existsFtoV(face))
  //     {
  //       ConstFaceVertices fv = topology->getFtoV(face);
  //       coord.reserve(fv.size());
  //       for (unsigned i=0; i<fv.size(); i++)
  //         coord.push_back(geometry->getPoint(fv[i]));
  //     }
  //   else if (topology->existsFtoE(face))
  //     {
  //       ConstFaceEdges fe = topology->getFtoE(face);
  //       for (unsigned i=0; i<fe.size(); i++)
  //         if (!topology->existsEtoV(fe[i]))
  //           throw Exception("faceCoord:  insufficient information");
  //       if (fe.size() == 0) // empty -- nothing to do
  //         ;
  //       else if (fe.size() == 1) // only one edge -- just get its coord
  //         {
  //           ConstEdgeVertices e1 = topology->getEtoV(fe[0]);
  //           coord.reserve(2);
  //           coord.push_back(geometry->getPoint(e1[0]));
  //           coord.push_back(geometry->getPoint(e1[1]));
  //         }
  //       else // Regular part
  //         {
  //           Id v1, v2, vn;
  //           ConstEdgeVertices e1 = topology->getEtoV(fe[0]);
  //           ConstEdgeVertices e2 = topology->getEtoV(fe[1]);
//#if 0
  //           Id e1_v0 =  e1[0];
  //           Id e1_v1 =  e1[1];
  //           Id e2_v0 =  e2[0];
  //           Id e2_v1 =  e2[1];
  //           if ( e1_v0.getValue() ) {}
  //           if ( e1_v1.getValue() ) {}
  //           if ( e2_v0.getValue() ) {}
  //           if ( e2_v1.getValue() ) {}
//#endif
  //           // Get v1 (the common vertex between e1 and e2)
  //           if (e1[0] == e2[0] || e1[0] == e2[1])
  //             v1 = e1[0];
  //           else if (e1[1] == e2[0] || e1[1] == e2[1])
  //             v1 = e1[1];
  //           else {
  //             std::cout << "e1[0] vert = " << e1[0] << ", coord " << geometry->getPoint(e1[0]) << std::endl;
  //             std::cout << "e1[1] vert = " << e1[1] << ", coord " << geometry->getPoint(e1[1]) << std::endl;
  //             std::cout << "e2[0] vert = " << e2[0] << ", coord " << geometry->getPoint(e2[0]) << std::endl;
  //             std::cout << "e2[1] vert = " << e2[1] << ", coord " << geometry->getPoint(e2[1]) << std::endl;
  //             Assert(false, "faceCoord:  incorrect information");
  //             //throw Exception("faceCoord:  incorrect information");
  //           }
  //           // Get the other two vertices of these edges
  //           v2 = Id(e2[0].toUnsigned() + e2[1].toUnsigned() - v1.toUnsigned()); // second vertex of face
  //           vn = Id(e1[0].toUnsigned() + e1[1].toUnsigned() - v1.toUnsigned()); // last vertex of face
  //           // Fill the coord vector by traversing through face edges
  //           coord.reserve(fe.size());
  //           coord.push_back(geometry->getPoint(v1));
  //           coord.push_back(geometry->getPoint(v2));
  //           for (unsigned i=2; i<fe.size()-1; i++)
  //             {
  //               ConstEdgeVertices e = topology->getEtoV(fe[i]);
  //               v2 = Id(e[0].toUnsigned() + e[1].toUnsigned() - v2.toUnsigned());
  //               coord.push_back(geometry->getPoint(v2));
  //             }
  //           coord.push_back(geometry->getPoint(vn));
  //         }
  //     }
  //   else
  //     throw Exception("faceCoord:  insufficient information");

  //    return;
  // }


  // Face Area

  // Vector3d faceArea(GeometryData *geometry, TopologyData *topology,
  //                     Id const &face)
  // {
  //   std::vector<Vector3d > coord;

  //    faceCoordinates(geometry, topology, face, coord);
  //   return faceArea(coord);
  // }




  // Vector3d faceArea(GeometryData *geometry, std::vector<Id> &vertices)
  // {
  //   std::vector<Vector3d > coord;
  //   coord.reserve(vertices.size());

  //    for (unsigned i=0; i<vertices.size(); i++)
  //     coord.push_back(geometry->getPoint(vertices[i]));

  //    return faceArea(coord);
  // }




  // Vector3d faceArea(GeometryData const * geometry, FaceVertices const &fv)
  // {
  //   std::vector<Vector3d > coord;

  //    coord.reserve(fv.size());

  //    for (unsigned i=0; i<fv.size(); i++)
  //     coord.push_back(geometry->getPoint(fv[i]));

  //    return faceArea(coord);
  // }




  Vector3d faceArea(std::vector<Vector3d > const &coord)
  {
    // Some of you may be wondering why I changed this routine.
    //   1. It's significantly faster than the old implementation.
    //   2. It gives exactly the same results. (Wayne Oaks)
    //
    // I am changing this back to how it was in 1.49.  Statement 2 above
    // is false.  Any faces which are offset significantly from the 
    // origin will suffer from severe floating point accuracy issues with
    // the 1.50 algorithm because it is not subtracting a suitable hub point
    // from the coordinates.


    // Sanity check
    if (coord.size() < 3)
      throw Exception("faceArea:  less than 3 vertices");
    else if (coord.size() == 3)
      return triangleArea(coord[0], coord[1], coord[2]);

    unsigned i;
    unsigned v1, v2;
    Vector3d fhp;   // face hub point
    Vector3d fArea; // total face area

    // Compute face hub point (vertex coordinate average);

    for (i=0; i<coord.size(); i++)
      fhp += coord[i];
    fhp *= (1.0/(double) coord.size()); // valid since size is non-zero

    // Loop over edges of face and accumulate area

    for (i=0; i<coord.size(); i++)
      {
        v1 = i;
        v2 = v1+1;
        if (v2 == coord.size())
	  v2 = 0;
        fArea += triangleArea(coord[v2],fhp,coord[v1]);
      }

    return fArea;
  }





  // Face Centroid

  // Vector3d faceCentroid(GeometryData *geometry, TopologyData *topology,
  //                         Id const &face, double areaTol)
  // {
  //   std::vector<Vector3d > coord;

  //    faceCoordinates(geometry, topology, face, coord);
  //   return faceCentroid(coord, areaTol);
  // }




  // Vector3d faceCentroid(GeometryData *geometry, std::vector<Id> &vertices,
  //                         double areaTol)
  // {
  //   std::vector<Vector3d > coord;

  //    coord.reserve(vertices.size());

  //    for (unsigned i=0; i<vertices.size(); i++)
  //     coord.push_back(geometry->getPoint(vertices[i]));

  //    return faceCentroid(coord, areaTol);
  // }




  // Vector3d faceCentroid(GeometryData *geometry, FaceVertices const &fv,
  //                         double areaTol)
  // {
  //   std::vector<Vector3d > coord;

  //    coord.reserve(fv.size());

  //    for (unsigned i=0; i<fv.size(); i++)
  //     coord.push_back(geometry->getPoint(fv[i]));

  //    return faceCentroid(coord, areaTol);
  // }




  Vector3d faceCentroid(std::vector<Vector3d > const &coord, double areaTol)
  {
    // Sanity check
    if (coord.size() < 3)
      throw Exception("faceCentroid:  less than 3 vertices");
    else if (coord.size() == 3)
      return triangleCentroid(coord[0], coord[1], coord[2]);

    unsigned i;
    Vector3d fhp; // Face hub point

    for (i=0; i<coord.size(); i++)
      fhp += coord[i];
    fhp *= (1.0/(double) coord.size()); // okay because size non-zero

    unsigned v1, v2;
    Vector3d centroid;
    Vector3d triArea;
    double triAreaMag;
    double areaMagTot = 0.0;

    // Loop over edges and accumulate area weighted centroid

    for (i=0; i<coord.size(); i++)
      {
        v1 = i;
        v2 = v1+1;
        if (v2 == coord.size())
	  v2 = 0;
        triArea = triangleArea(coord[v2], fhp, coord[v1]);
        triAreaMag = triArea.length();
        centroid += triangleCentroid(coord[v2],fhp,coord[v1])*triAreaMag;
        areaMagTot += triAreaMag;
      }

    if (areaMagTot < areaTol) // 0 area -- return vertex average!
      {
        centroid = 0.0;
        for (i=0; i<coord.size(); i++)
          centroid += coord[i];
        centroid *= (1.0/(double) coord.size());
        return centroid;
      }

    centroid *= (1.0/areaMagTot);
    return centroid;
  }




  // Face area and Centroid

  // void faceAreaCentroid(GeometryData *geometry, TopologyData *topology,
  //                       Id const &face,
  //                       Vector3d &area, Vector3d &centroid, double areaTol)
  // {
  //   std::vector<Vector3d > coord;

  //    faceCoordinates(geometry, topology, face, coord);
  //   faceAreaCentroid(coord, area, centroid, areaTol);

  //    return;
  // }




  // void faceAreaCentroid(GeometryData *geometry, std::vector<Id> &vertices,
  //                       Vector3d &area, Vector3d &centroid, double areaTol)
  // {
  //   std::vector<Vector3d > coord;

  //    coord.reserve(vertices.size());

  //    for (unsigned i=0; i<vertices.size(); i++)
  //     coord.push_back(geometry->getPoint(vertices[i]));

  //    faceAreaCentroid(coord, area, centroid, areaTol);

  //    return;
  // }




  // void faceAreaCentroid(GeometryData *geometry, FaceVertices const &fv,
  //                       Vector3d &area, Vector3d &centroid, double areaTol)
  // {
  //   std::vector<Vector3d > coord;

  //    coord.reserve(fv.size());

  //    for (unsigned i=0; i<fv.size(); i++)
  //     coord.push_back(geometry->getPoint(fv[i]));

  //    faceAreaCentroid(coord, area, centroid, areaTol);

  //    return;
  // }




  void faceAreaCentroid(std::vector<Vector3d > const &coord,
                        Vector3d &area, Vector3d &centroid, double areaTol)
  {
    if (coord.size() < 3)
      throw Exception("faceAreaCentroid:  less than 3 vertices");
    else if (coord.size() == 3) {
      triangleAreaCentroid(coord[0], coord[1], coord[2], area, centroid);
      return;
    }
    
    unsigned i;
    Vector3d fhp; // Face hub point

    for (i=0; i<coord.size(); i++)
      fhp += coord[i];
    fhp *= (1.0/(double) coord.size()); // okay because size non-zero

    unsigned v1, v2;
    Vector3d triArea;
    double triAreaMag;
    double areaMagTot = 0.0;

    // Loop over edges and accumulate area and area weighted centroid

    area     = 0.0; // initialize for accumulation
    centroid = 0.0;

    for (i=0; i<coord.size(); i++)
      {
        v1 = i;
        v2 = v1+1;
        if (v2 == coord.size())
	  v2 = 0;
        triArea = triangleArea(coord[v2], fhp, coord[v1]);
        area += triArea;
        triAreaMag = triArea.length();
        centroid += triangleCentroid(coord[v2],fhp,coord[v1])*triAreaMag;
        areaMagTot += triAreaMag;
      }

    if (areaMagTot < areaTol) // 0 area -- return vertex average!
      {
        centroid = 0.0;
        for (i=0; i<coord.size(); i++)
          centroid += coord[i];
        centroid *= (1.0/(double) coord.size());
        return;
      }

    centroid *= (1.0/areaMagTot);

    return;
  }


  // double facePointDistanceSquared(GeometryData *geometry, 
  //                                 TopologyData *topology,
  //                                 Id const &face,
  //                                 Vector3d const &querypt)
  // {
  //   std::vector<Vector3d > coord;
  //       
  //   faceCoordinates(geometry, topology, face, coord);
  //   return facePointDistanceSquared(coord, querypt);
  // }


  // double facePointDistanceSquared(GeometryData *geometry, 
  //                                 FaceVertices const &fv,
  //                                 Vector3d const &querypt)
  // {
  //   if (fv.size() == 3)
  //     return trianglePointDistanceSquared(geometry, fv, querypt);

  //    std::vector<Vector3d > coord;
  //       
  //   coord.reserve(fv.size());
  //       
  //   for (unsigned i=0; i<fv.size(); i++)
  //     coord.push_back(geometry->getPoint(fv[i]));
  //       
  //   return facePointDistanceSquared(coord, querypt);        
  // }


  double facePointDistanceSquared(std::vector<Vector3d > const &coord,
                                  Vector3d const &querypt)
  {
    if (coord.size() < 3)
      throw Exception("facePointDistanceSquared: less than 3 vertices");
    else if (coord.size() == 3)
      return trianglePointDistanceSquared(coord[0], coord[1], coord[2],
                                          querypt);
        
        
    // find minimum distance to any triangle formed by an edge and
    // hub point
    // fixme:  THIS MAY NOT WORK FOR CONCAVE FACES!
    double d_min = std::numeric_limits<double>::max();
    Vector3d centroid(faceCentroid(coord));
    unsigned i, v1, v2;
    for (i=0; i<coord.size(); i++)
      {
        v1 = i;
        v2 = v1+1;
        if (v2 == coord.size())
	  v2 = 0;
        d_min = 
          std::min(d_min, trianglePointDistanceSquared(coord[v2], 
                                                       centroid,
                                                       coord[v1], 
                                                       querypt));
      }
    return d_min;
  }






  // Face Planarity

  // double facePlanarity(GeometryData *geometry, TopologyData *topology,
  //        Id const &face)
  // {
  //   std::vector<Vector3d > coord;
  //   
  //   faceCoordinates(geometry, topology, face, coord);
  //   return facePlanarity(coord);
  // }
  
  
  
  
  // double facePlanarity(GeometryData *geometry, std::vector<Id> &vertices)
  // {
  //   std::vector<Vector3d > coord;

//    coord.reserve(vertices.size());

//    for (unsigned i=0; i<vertices.size(); i++)
//     coord.push_back(geometry->getPoint(vertices[i]));

//    return facePlanarity(coord);
// }




// double facePlanarity(GeometryData *geometry, FaceVertices const &fv)
// {
//   std::vector<Vector3d > coord;

//    coord.reserve(fv.size());

//    for (unsigned i=0; i<fv.size(); i++)
//     coord.push_back(geometry->getPoint(fv[i]));

//    return facePlanarity(coord);
// }




  double facePlanarity(std::vector<Vector3d > const &coord)
  {
    // Sanity check
    if (coord.size() < 3)
      throw Exception("facePlanarity:  less than 3 vertices");
    else if (coord.size() == 3)
      return 1.0;

    //get face area vector;
    Vector3d faceArea,faceCentroid;
    try {
      faceAreaCentroid(coord,faceArea,faceCentroid);
    }
    catch (...) {
      return 0.0;
    }

    unsigned i,ip;
    double a_i;
    double a_tot=0.0;
    Vector3d u,v;
    for (i=0; i<coord.size(); i++) {
      ip=i+1; if(ip==coord.size()) ip=0;
      u=coord[i ]-faceCentroid;
      v=coord[ip]-faceCentroid;
      a_i = ((u.cross(v)).mag());
      a_tot+=a_i;
    }
    a_tot *=0.5;

    //check
    if (a_tot==0.0) return 0.0;

    return faceArea.mag()/a_tot;

  }

  double facePlanarity(std::vector<Vector3d > const &coord,
		       Vector3d const &faceArea,
		       Vector3d const &faceCentroid)
  {
    // Sanity check
    if (coord.size() < 3)
      throw Exception("facePlanarity:  less than 3 vertices");
    else if (coord.size() == 3)
      return 1.0;

    unsigned i,ip;
    double a_i;
    double a_tot=0.0;
    Vector3d u,v;
    for (i=0; i<coord.size(); i++) {
      ip=i+1; if(ip==coord.size()) ip=0;
      u=coord[i ]-faceCentroid;
      v=coord[ip]-faceCentroid;
      a_i = ((u.cross(v)).mag());
      a_tot+=a_i;
    }
    a_tot *=0.5;

    //check
    if (a_tot==0.0) return 0.0;

    return faceArea.mag()/a_tot;

  }



  // Face twist

  // double faceTwist(GeometryData *geometry, TopologyData *topology,
  //        Id const &face)
  // {
  //   std::vector<Vector3d > coord;
  //   
  //   faceCoordinates(geometry, topology, face, coord);
  //   return faceTwist(coord);
  // }
  
  
  
  
  // double faceTwist(GeometryData *geometry, std::vector<Id> &vertices)
  // {
  //   std::vector<Vector3d > coord;

//    coord.reserve(vertices.size());

//    for (unsigned i=0; i<vertices.size(); i++)
//     coord.push_back(geometry->getPoint(vertices[i]));

//    return faceTwist(coord);
// }




// double faceTwist(GeometryData *geometry, FaceVertices const &fv)
// {
//   std::vector<Vector3d > coord;

//    coord.reserve(fv.size());

//    for (unsigned i=0; i<fv.size(); i++)
//     coord.push_back(geometry->getPoint(fv[i]));

//    return faceTwist(coord);
// }




  double faceTwist(std::vector<Vector3d > const &coord)
  {
    // Sanity check
    if (coord.size() < 3)
      throw Exception("faceTwist:  less than 3 vertices");
    else if (coord.size() == 3)
      return 0.0;

    //get face area vector;
    Vector3d faceArea, faceCentroid;
    try {
      faceAreaCentroid(coord,faceArea,faceCentroid);
      faceArea.normalize();
    }
    catch (...) {
      return 1.0;
    }

    return faceTwist(coord,faceArea,faceCentroid);
  }

  double faceTwist(std::vector<Vector3d > const &coord,
		   Vector3d const &faceArea,
		   Vector3d const &faceCentroid)
  {
    // Sanity check
    if (coord.size() < 3)
      throw Exception("faceTwist:  less than 3 vertices");
    else if (coord.size() == 3)
      return 0.0;


    unsigned i,ip;
    double a;
    double a_tot=0.0;
    double twist = 0.0; 
    std::vector< Vector3d > T;
    T.resize(coord.size());
    Vector3d u,v;

    for (i=0; i<coord.size(); i++) {
      ip=i+1; if(ip==coord.size()) ip=0;
      u=coord[i ]-faceCentroid;
      v=coord[ip]-faceCentroid;
      T[i] = u.cross(v);
    }

    
    for (i=0; i<coord.size(); i++) {

      a = T[i].dot(faceArea);
      
      if (a<0.0) {
	twist -= a;
	a_tot -= a;
      }
      else {
	a_tot += a;
      }
    }

    //check
    if (a_tot==0.0) return 1.0;
    
    //total twist
    return twist/a_tot;
    
  }
}
