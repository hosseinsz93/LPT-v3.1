// This file and its contents are the property of adapco, Ltd. and contain
// confidential and proprietary information.  The unauthorized use,
// distribution, or duplication of this information is prohibited.
// All rights reserved.
//
// Copyright (C) 2004 by adapco, Ltd.

#include "GeomEdgeCalculations.h"

#include "GeomVector.h"

// #include "TopologyData.h"
// #include "GeometryData.h"

namespace Geom
{
//   void edgeVertices(TopologyData *topology, Id const &edge,
//                     Id &v1, Id &v2)
//   {
//     if (topology->existsEtoV(edge))
//       {
//         ConstEdgeVertices ev = topology->getEtoV(edge);
//         v1 = ev[0];
//         v2 = ev[1];
//       }
//     else
//       throw Exception("edgeVertices:  insufficient information");
//   }
  // get the endpoint coordinates of an edge
//   void edgeCoordinates(GeometryData *geometry, TopologyData *topology, 
//                        Id const &edge, Vector3d &coord1, 
//                        Vector3d &coord2)
//   {
//     Id v1, v2;
//     edgeVertices(topology, edge, v1, v2);
//     coord1 = geometry->at(v1);
//     coord2 = geometry->at(v2);
//   }

//   double edgePointDistanceSquared(GeometryData *geometry,
//                                   TopologyData *topology,
//                                   Id const &edge,
//                                   Vector3d const &querypt)
//   {
//     Vector3d coord1(0), coord2(0);
//     edgeCoordinates(geometry, topology, edge, coord1, coord2);
//     return edgePointDistanceSquared(coord1, coord2, querypt);
//   }
        
//   double edgePointDistanceSquared(GeometryData *geometry,
//                                   ConstEdgeVertices const &ev,
//                                   Vector3d const &querypt)
//   {
//     if (ev.size() != 2)
//       throw Exception("edgePointDistanceSquared: not a simple edge");
//     return edgePointDistanceSquared(geometry->at(ev[0]),
//                                     geometry->at(ev[1]),
//                                     querypt);       
//   }

//   double edgePointDistanceSquared(GeometryData *geometry,
//                                   Id const &v1, Id const &v2,
//                                   Vector3d const &querypt)
//   {
//     return edgePointDistanceSquared(geometry->at(v1),
//                                     geometry->at(v2),
//                                     querypt);
//   }

  double edgePointDistanceSquared(Vector3d const &coord1,
                                  Vector3d const &coord2,
                                  Vector3d const &querypt,
				  Vector3d *closestPoint)
  {
    // distance from point to segment
    Vector3d diff(querypt - coord1);
    Vector3d dir(coord2 - coord1);
    double t = diff.dot(dir);

    if (t > 0.0)
      {
        double len_squared = dir.mag2();
        if (t >= len_squared)
	{
          diff -= dir;
	  if (closestPoint)
	      *closestPoint = coord2;
	}
        else
	{
            t /= len_squared;
            diff -= t*dir;
	    if (closestPoint)
		*closestPoint = querypt - diff;
	}
      }
    else if (closestPoint)
	*closestPoint = coord1;
    return diff.mag2();
  }



  Vector3d edgeProjectPoint(Vector3d const &coord1,
                            Vector3d const &coord2,
                            Vector3d const &coord)
  {
    // Project coord onto edge coord1-coord2

    Vector3d vec1(coord2 - coord1);
    Vector3d vec2(coord  - coord1);

    double dot = vec1.dot(vec2);
    if (dot > 0.0)
      {
        double lengthsq = vec1.mag2();
        if (dot >= lengthsq)
          {
            return coord2;
          }
        else
          {
            double ratio = dot/lengthsq;
            return Vector3d(coord1 + ratio*vec1);
          }
      }

    return coord1;
  }



  double edgeProjectRatio(Vector3d const &coord1,
                          Vector3d const &coord2,
                          Vector3d const &coord)
  {
    Vector3d vec1(coord2 - coord1);
    Vector3d vec2(coord  - coord1);

    double dot = vec1.dot(vec2);
    if (dot <= 0.0)
      return 0.0;

    double lengthsq = vec1.mag2();
    if (dot >= lengthsq)
      return 1.0;

    return dot/lengthsq;
  }



  double edgeProjectRatioExtended(Vector3d const &coord1,
                                  Vector3d const &coord2,
                                  Vector3d const &coord)
  {
    Vector3d vec1(coord2 - coord1);
    Vector3d vec2(coord  - coord1);
    double dot = vec1.dot(vec2);
    double lengthsq = vec1.mag2();
    return dot/lengthsq;
  }


  
  Vector3d edgeCenter(Vector3d const &A, Vector3d const &B) 
  {
    return Vector3d((A + B) * 0.5);
  } 

  Vector3d edgeDirectional(Vector3d const &A, Vector3d const &B) 
  {  
    Vector3d r(B - A); 
    r.normalize();
    return r;
  }

}
