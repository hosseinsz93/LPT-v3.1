#include "GeomXFaceCalculations.h"

#include <cassert>

#include "GeomBoundingBox.h"
#include "GeomXCoordSystem.h"
#include "GeomXUtilities.h"

// #include "GeomTopologyDataPointers.h"
// #include "GeomXVertexCoordInterface.h"

#include "GeomXEdgeCalculations.h"
#include "GeomXException.h"
// #include "GeomXGenericVertices_t.h"


namespace GeomX
{
#if 0
  void faceCoordinates(MKXCore::VertexCoordInterface const * geometry,
                       MKCore::TopologyData const * topology,
                       MKCore::Id const &face,
                       std::vector< Vector3d > &coords)
  {
    coords.clear(); // QUESTION: Is this expensive or simply setting size=0?

    if (topology->existsFtoV(face))
      {
        MKCore::ConstFaceVertices fv = topology->getFtoV(face);
        coords.reserve(fv.size());
        for (unsigned i=0; i<fv.size(); i++)
          coords.push_back((*geometry)[fv[i]]);
      }
    else if (topology->existsFtoE(face))
      {
        MKCore::ConstFaceEdges fe = topology->getFtoE(face);
        for (unsigned i=0; i<fe.size(); i++)
          if (!topology->existsEtoV(fe[i]))
            throw Exception("faceCoordinates: insufficient information");
        if (fe.size() == 0) // empty -- nothing to do
          ;
        else if (fe.size() == 1) // only one edge -- just get its coord
          {
            MKCore::ConstEdgeVertices e1 = topology->getEtoV(fe[0]);
            coords.reserve(2);
            coords.push_back((*geometry)[e1[0]]);
            coords.push_back((*geometry)[e1[1]]);
          }
        else // Regular part
          {
            MKCore::Id v1, v2, vn;
            MKCore::ConstEdgeVertices e1 = topology->getEtoV(fe[0]);
            MKCore::ConstEdgeVertices e2 = topology->getEtoV(fe[1]);
            // Get v1 (the common vertex between e1 and e2)
            if (e1[0] == e2[0] || e1[0] == e2[1])
              v1 = e1[0];
            else if (e1[1] == e2[0] || e1[1] == e2[1])
              v1 = e1[1];
            else
              throw Exception("faceCoordinates: incorrect information");
            // Get the other two vertices of these edges
            v2 = MKCore::Id(e2[0].toUnsigned() + e2[1].toUnsigned() - v1.toUnsigned()); // second vertex of face
            vn = MKCore::Id(e1[0].toUnsigned() + e1[1].toUnsigned() - v1.toUnsigned()); // last vertex of face
            // Fill the coord vector by traversing through face edges
            coords.reserve(fe.size());
            coords.push_back((*geometry)[v1]);
            coords.push_back((*geometry)[v2]);
            for (unsigned i=2; i<fe.size()-1; i++)
              {
                MKCore::ConstEdgeVertices e = topology->getEtoV(fe[i]);
                v2 = MKCore::Id(e[0].toUnsigned() + e[1].toUnsigned() - v2.toUnsigned());
                coords.push_back((*geometry)[v2]);
              }
            coords.push_back((*geometry)[vn]);
          }
      }
    else
      throw Exception("faceCoordinates: insufficient information");

    return;
  }
#endif

#if 0
  inline void helper(MKXCore::VertexCoordInterface const * geometry,
                     std::vector< Vector3d > &coords,
                     MKCore::Id const & id, MKCore::Id const & vId,
                     Vector3d const & c)
  {
    if (id == vId)
      coords.push_back(c);
    else
      coords.push_back((*geometry)[id]);
  }
#endif

#if 0
  void getFaceCoords(MKXCore::VertexCoordInterface const * geom,
                     MKCore::FaceVertices const & fv,
                     std::vector < Vector3d > & coords)
  {
    coords.resize(0);  // Clear without memory delete
    coords.reserve(fv.size());
    for (unsigned i=0; i<fv.size(); i++)
      coords.push_back((*geom)[fv[i]]);
  }
#endif

  
#if 0
  Vector3d faceArea(MKXCore::VertexCoordInterface const * geom,
                    MKCore::FaceVertices const & fv)
  {
    std::vector < Vector3d > coord;
    coord.reserve(fv.size());
    for(unsigned i=0; i<fv.size(); i++)
      coord.push_back((*geom)[fv[i]]);
    return MKCore::faceArea(coord);
  }
#endif


  Vector3d triangleAreaCoords(Vector3d const & p,
                              Vector3d const & c0,
                              Vector3d const & c1,
                              Vector3d const & c2)
  {
    Vector3d c[3];
    c[0] = c0;
    c[1] = c1;
    c[2] = c2;
    return triangleAreaCoords(p, c);
  }


  Vector3d triangleAreaCoords(Vector3d const & p,
                              Vector3d * cs)
  {
    Vector3d areaCoord;
    Vector3d triArea(Geom::triangleArea(cs[0], cs[1], cs[2]));
    double area = triArea.length();

    for (int i=0; i<3; ++i)
    {
      Vector3d tmp(cs[i]);
      cs[i] = p;
      Vector3d partArea(Geom::triangleArea(cs[0], cs[1], cs[2]));
      double pArea = partArea.length();
      areaCoord(i) = sign(pArea/area, triArea.dot(partArea));
      cs[i] = tmp;
    }
    return areaCoord;
  }


  Vector3d triangleAreaCoords2D(Vector3d const & p,
                                Vector3d * cs)
  {
    double area = Geom::triangleArea2D(cs[0], cs[1], cs[2]);

    Vector3d areaCoord;
    for (int i=0; i<3; ++i)
    {
      Vector3d tmp(cs[i]);
      cs[i] = p;
      double pArea = Geom::triangleArea2D(cs[0], cs[1], cs[2]);
      areaCoord(i) = sign(pArea/area, area*pArea);
      cs[i] = tmp;
    }
    return areaCoord;
  }


  void triangleAreaCoordsConstrain(Vector3d const & p,
                                   Vector3d const & c0,
                                   Vector3d const & c1,
                                   Vector3d const & c2,
                                   Vector3d       & areaCoord)
  {
    Vector3d c[3];
    c[0] = c0;
    c[1] = c1;
    c[2] = c2;
    return triangleAreaCoordsConstrain(p, c, areaCoord);
  }


  void triangleAreaCoordsConstrain(Vector3d const & p,
                                   Vector3d c[],
                                   Vector3d & areaCoord)
  {
    unsigned u, v;
    int negCnt=0, negIdx=-1, nonNegIdx=-1;

    for (u=0; u<3; ++u)
    {
      if (areaCoord(u) < 0.0)
      {
        ++negCnt;
        negIdx = u;
        areaCoord(u) = 0.0;
      }
      else
      {
        nonNegIdx = u;
        if (areaCoord(u) > 1.0)
          areaCoord(u) = 1.0;
      }
    }


    if (negCnt == 0)
      return;

    if (negCnt == 2)
    {
      areaCoord(nonNegIdx) = 1.0;
      return;
    }

    assert(negCnt == 1);
    u = negIdx + 1;
    if (u >= 3) u = 0;
    v = u + 1;
    if (v >= 3) v = 0;

    areaCoord(v) = lineProjectPointPnt_t(p, c[u], c[v]);
    if (areaCoord(v) < 0.0) areaCoord(v) = 0.0;
    if (areaCoord(v) > 1.0) areaCoord(v) = 1.0;
    areaCoord(u) = 1.0 - areaCoord(v);
  }


  void makeCoordsTrans(std::vector < Vector3d > const & coords,
                       CoordSystem & cs)
  {
    Vector3d area, centroid;
    Geom::faceAreaCentroid(coords, area, centroid);
    cs.setOrthoNormalSet(centroid, area);
  }


  void transformCoords(CoordSystem & cs,
                       std::vector < Vector3d > const & coords,
                       std::vector < Vector3d > & coords_cs)
  {

    std::vector < Vector3d >().swap(coords_cs);
    coords_cs.reserve(coords.size());
    std::vector < Vector3d >::const_iterator i;
    for (i=coords.begin(); i!=coords.end(); ++i)
      coords_cs.push_back(cs.transInto(*i));
  }


  // All the coordinates must be in the coordinate system defined
  // by coords.
  bool facePierced(std::vector < Vector3d > const & coords,
                   Vector3d const & e0, Vector3d const & e1)
  {
    // If on the same side of the plane, there is no overlap.
    double c = e0(2) * e1(2);
    if (c >= 0.0)
      return false;

    // Find the intersection point in the xy plane.
    Vector3d p;
    if (c == 0.0)
      p = (e0(2) == 0.0) ? e0 : e1;
    else
    {
      double ratioa = e1(2) / (e1(2) - e0(2));
      double ratiob = e0(2) / (e1(2) - e0(2));
      p = Vector3d(ratioa * e0(0) - ratiob * e1(0),
                   ratioa * e0(1) - ratiob * e1(1),
                   0.0);
    }
    return pointInsideOnFace(coords, p);
  }


  // Are the corner points closer to and inside this face.
  bool facePointDist(std::vector < Vector3d > const & coords,
                     Vector3d const & pnt, double & dist)
  {
    // If the distance from the plane is not less than the currently
    // found minimum distance, no more need to check.
    if (std::abs(pnt(2)) >= dist)
      return false;

    // If the point is inside the face, use its distance as the new minimum.
    if (pointInsideOnFace(coords, pnt))
      dist = std::abs(pnt(2));

    return (dist == 0.0);
  }


  double faceFaceDistance(std::vector < Vector3d > const & coordsA,
                          std::vector < Vector3d > const & coordsB,
                          double tol, double * maxFacesDistance)
  {
    // If the distance between bounding boxes is > maxFacesDistance,
    // return a very large distance.
    if (maxFacesDistance)
    {
      Geom::BoundingBox bb0; bb0.expand(coordsA.begin(), coordsA.end());
      Geom::BoundingBox bb1; bb1.expand(coordsB.begin(), coordsB.end());
      if (bb0.distance(bb1) > *maxFacesDistance)
        return std::numeric_limits<double>::max();
    }

    // check distances between edges.
    double result = std::numeric_limits<double>::max();

    CoordSystem csA;
    std::vector < Vector3d > coordsA_csA, coordsB_csA;
    makeCoordsTrans(coordsA, csA);
    transformCoords(csA, coordsA, coordsA_csA);
    transformCoords(csA, coordsB, coordsB_csA);

    CoordSystem csB;
    std::vector < Vector3d > coordsB_csB, coordsA_csB;
    makeCoordsTrans(coordsB, csB);
    transformCoords(csB, coordsB, coordsB_csB);
    transformCoords(csB, coordsA, coordsA_csB);


    int a0, a1, b0, b1;
    for (a0=0, a1=1; a0<(int)coordsA.size(); ++a0, ++a1)
    {
      if (a1 >= (int)coordsA.size()) a1 = 0;

      for (b0=0, b1=1; b0<(int)coordsB.size(); ++b0, ++b1)
      {
        if (b1 >= (int)coordsB.size()) b1 = 0;

        double dist = edgeEdgeDistancePnt_t(coordsA[a0], coordsA[a1],
                                            coordsB[b0], coordsB[b1], tol);
        result = std::min(result, dist);
        if (result < tol)
          return 0.0;


        // See if either edge pierces the other face.
        if (facePierced(coordsA_csA, coordsB_csA[b0], coordsB_csA[b1]))
          return 0.0;
        if (facePierced(coordsB_csB, coordsA_csB[a0], coordsA_csB[a1]))
          return 0.0;
      }
    }


    // Are the corner points closer to and inside the other face.
    for (a0=0; a0<(int)coordsA.size(); ++a0)
    {
      if(facePointDist(coordsB_csB, coordsA_csB[a0], result))
        return 0.0;
      if (result < tol)
        return 0.0;
    }

    for (b0=0; b0<(int)coordsB.size(); ++b0)
    {
      if(facePointDist(coordsA_csA, coordsB_csA[b0], result))
        return 0.0;
      if (result < tol)
        return 0.0;
    }


    if (maxFacesDistance)
    {
      return (result <= *maxFacesDistance) ? result
                                           : std::numeric_limits<double>::max();
    }

    return result;
  }


  bool faceIsSelfIntersecting(std::vector < Vector3d > const & cs)
  {
    int otn = cs.size();

    // Triangles and less can't intersect.
    if (otn <= 3)
      return false;

    int inn = otn - 3;
    for (int i=0; i<otn-2; ++i)
    {
      for (int j=0; j<inn; ++j)
      {
        // Process edge i with i+j+2
        int ks = i + j + 2;
        int ke = ks + 1;
        if (ke >= otn)
          ke = 0;

        bool result =
          edgeEdgeIntersect(toVector2D(cs[i]),  toVector2D(cs[i+1]),
                            toVector2D(cs[ks]), toVector2D(cs[ke]),
                            almostZero());
        if (result)
          return true;
      }
    }
    return false;
  }
} // namespace MKXGeom


// Automatic setting of emacs local variables.
// Local Variables:
// mode: C++
// tab-width: 8
// End:
