#include "GeomXEdgeCalculations.h"

#include <limits>

// #include "common/Tensor.h"
// #include "common/TensorAlgorithms.h"

namespace GeomX
{
  template < typename T >
    inline T square(T const &a) { return a*a; }


#if 0
  // get the endpoint coordinates of an edge
  void edgeCoordinates(VertexCoordInterface * geometry,
                       Geom::TopologyData *topology, Geom::Id const &edge,
                       Vector3d &coord1, Vector3d &coord2)
  {
    MKCore::Id v1, v2;
    MKCore::edgeVertices(topology, edge, v1, v2);
    coord1 = (*geometry)[v1];
    coord2 = (*geometry)[v2];
  }
#endif
  

  bool edgeEdgeIntersect(Vector2d const &p0,
                         Vector2d const &p1,
                         Vector2d const &q0,
                         Vector2d const &q1,
                         double &alphap, double &alphaq)
  {
    // Normal to q0, q1.
    Vector2d nq(q0(1)-q1(1), q1(0)-q0(0));
    // See if p line segment crosses q0->q1.
    double wec_p0 = Vector2d(p0-q0).dot(nq);
    double wec_p1 = Vector2d(p1-q0).dot(nq);
    if (wec_p0*wec_p1 <= 0.0)
    {
      // Normal to q0, q1.
      Vector2d np(p0(1)-p1(1), p1(0)-p0(0));
      // See if q line segment crosses p0->p1.
      double wec_q0 = Vector2d(q0-p0).dot(np);
      double wec_q1 = Vector2d(q1-p0).dot(np);
      if (wec_q0*wec_q1 <= 0.0)
      {
        alphap = wec_p0 / (wec_p0 - wec_p1);
        alphaq = wec_q0 / (wec_q0 - wec_q1);
        return true;
      }
    }

    return false;
  }


  bool edgeEdgeIntersect(Vector2d p0, Vector2d p1,
                         Vector2d q0, Vector2d q1,
                         double linearTol)
  {
    Vector2d xAxis(p1 - p0);
    double pLeng = xAxis.normalize();

    Vector2d nvq(q1 - q0);
    double qLeng = nvq.normalize();

    if (pLeng < qLeng)
    {
      std::swap(p0, q0);
      std::swap(p1, q1);
      std::swap(xAxis, nvq);
      std::swap(pLeng, qLeng);
    }

    if (pLeng <= linearTol)
      return false;

    Vector2d yAxis(-xAxis(1), xAxis(0));

    q0 -= p0;
    q0 = Vector2d(xAxis.dot(q0), yAxis.dot(q0));
    if (std::abs(q0(1)) <= linearTol &&
        q0(0) >= -linearTol &&
        q0(0) <= pLeng + linearTol)
      return true;

    q1 -= p0;
    q1 = Vector2d(xAxis.dot(q1), yAxis.dot(q1));
    if (std::abs(q1(1)) <= linearTol &&
        q1(0) >= -linearTol &&
        q1(0) <= pLeng + linearTol)
      return true;

    if (q0(1) * q1(1) > 0.0)
      return false;

    double x = (q0(0) * q1(1) - q1(0) * q0(1)) / (q1(1) - q0(1));
    if (x < -linearTol)
      return false;
    if (x > pLeng + linearTol)
      return false;

    return true;
  }


#if 0
  Vector3d edgeEdgeIntersect(Vector3d p0, Vector3d p1,
                             Vector3d q0, Vector3d q1)
  {
    Vector3d pv(p1-p0);
    Vector3d qv(q0-q1);
    Tensor<3,double> A(Geom::Vector<3,Vector3d >(pv, qv, pv.cross(qv)));
    A = A.transpose();
    Vector3d b(q0-p0);
    Vector3d s(TensorAlgorithms::solvePVT(A, b));
    return Vector3d(0.5 * (p0 + s(0)*pv + q0 - s(1)*qv));
  }
#endif
  



  // Routines after this comment to main were added around 2009-Feb-25

  template < int N, typename T >
  double linePointDistance(VectorNT const &coord1,
                           VectorNT const &coord2,
                           VectorNT const &querypt)
  {
    // distance from point to line
    VectorNT diff(querypt - coord1);
    VectorNT dir(coord2 - coord1);
    double t = diff.dot(dir);
    double len_squared = dir.mag2();
    t /= len_squared;
    diff -= t*dir;
    return diff.length();
  }
  template double linePointDistance(Vector2d const &,
                                    Vector2d const &,
                                    Vector2d const &);
  template double linePointDistance(Vector3d const &,
                                    Vector3d const &,
                                    Vector3d const &);


  template < int N, typename T >
  T lineProjectPointVec_t(VectorNT const & pnt,
                          VectorNT const & p0, VectorNT const & pv)
  {
    T leng2 = pv.mag2();
    if (leng2 == 0.0)
      return 0.0;

    T ratio = pv.dot(VectorNT(pnt-p0)) / leng2;
    return ratio;
  }
  template double lineProjectPointVec_t(Vector2d const &,
                                        Vector2d const &,
                                        Vector2d const &);
  template double lineProjectPointVec_t(Vector3d const &,
                                        Vector3d const &,
                                        Vector3d const &);


  template < int N, typename T >
  inline VectorNT unit(VectorNT const & v)
  {
    T const l = v.length();
    return (l > T(0)) ? VectorNT(v/l) : v;
  }


  // The data comes in as  |  |
  //                       p  q
  //                       |  |
  template < typename T >
  T determinant(Geom::Vector<2,T> const & p, Geom::Vector<2,T> const & q)
  {
    T det = p(0)*q(1) - q(0)*p(1);
    return det;
  }


  // The data comes in as  |  |  |
  //                       p  q  n
  //                       |  |  |
  template < typename T >
  T determinant(Geom::Vector<3,T> const & p,
                Geom::Vector<3,T> const & q, Geom::Vector<3,T> const & r)
  {
    // p(0) * | q(1) r(1) |
    //        | q(2) r(2) |

    // p(1) * | q(2) r(2) |
    //        | q(0) r(0) |

    // p(2) * | q(0) r(0) |
    //        | q(1) r(1) |

    T det = p(0) * determinant(Geom::Vector<2,T>(q(1),q(2)), Geom::Vector<2,T>(r(1),r(2)))
          + p(1) * determinant(Geom::Vector<2,T>(q(2),q(0)), Geom::Vector<2,T>(r(2),r(0)))
          + p(2) * determinant(Geom::Vector<2,T>(q(0),q(1)), Geom::Vector<2,T>(r(0),r(1)));
    return det;
  }


  // Calculates the distance between two lines.  This routine does not handle
  // segments that do not cross.
  template < int N, typename T >
  T lineLineDistance(VectorNT const & p0, VectorNT const & p1,
                     VectorNT const & q0, VectorNT const & q1)
  {
    Geom::Vector<3,T> pv(p1 - p0);
    Geom::Vector<3,T> qv(q1 - q0);
    Geom::Vector<3,T> n(pv.cross(qv));
    n.normalize();
    T det = determinant(pv, qv, n);
    Geom::Vector<3,T> pq(q0 - p0);
    T dist = determinant(pv, qv, pq) / det;
    return std::abs(dist);
  }


  // Line is used in the title since they are assumed to be infinitely long.
  template < int N, typename T >
  void lineLineParallelPoint(VectorNT const &p0, VectorNT const &p1,
                             VectorNT const &q0, VectorNT const &q1,
                             VectorNT const &pv, VectorNT const &qv,
                             double &ratiop, double &ratioq)
  {
    if (pv.length() >= qv.length())
    {
      T ratio0 = lineProjectPointVec_t(q0, p0, pv);
      T ratio1 = lineProjectPointVec_t(q1, p0, pv);
      if (std::abs(ratio0) < std::abs(ratio1))
      {
        ratiop = ratio0;
        ratioq = 0.0;
      }
      else
      {
        ratiop = ratio1;
        ratioq = 1.0;
      }
    }
    else
    {
      T ratio0 = lineProjectPointVec_t(p0, q0, qv);
      T ratio1 = lineProjectPointVec_t(p1, q0, qv);
      if (std::abs(ratio0) < std::abs(ratio1))
      {
        ratioq = ratio0;
        ratiop = 0.0;
      }
      else
      {
        ratioq = ratio1;
        ratiop = 1.0;
      }
    }
  }


  template < typename T >
  void lineLineIntersectPnt_t(Geom::Vector<2,T> const &p0, Geom::Vector<2,T> const &p1,
                              Geom::Vector<2,T> const &q0, Geom::Vector<2,T> const &q1,
                              T & ratiop, T & ratioq)
  {
    Geom::Vector<2,T> pv(p1 - p0);
    Geom::Vector<2,T> qv(q1 - q0);
    // For very small angles, set det to 0 to force the parallel intersector.
    T sin = std::abs(determinant(unit(pv), unit(qv)));
    T det = (sin < 0.01) ? 0.0 : determinant(pv, qv);

    if (det != 0.0)
    {
      Geom::Vector<2,T> pq(q0 - p0);
      ratiop = determinant(pq, qv) / det;
      ratioq = -determinant(pv, pq) / det;
      return;
    }

    lineLineParallelPoint(p0, p1, q0, q1, pv, qv, ratiop, ratioq);
  }

  // Instantiate
  template void lineLineIntersectPnt_t(Vector2d const &p0, Vector2d const &p1,
                                       Vector2d const &q0, Vector2d const &q1,
                                       double & ratiop, double & ratioq);


  template < typename T >
  void lineLineIntersectPnt_t(Geom::Vector<3,T> const & p0, Geom::Vector<3,T> const & p1,
                              Geom::Vector<3,T> const & q0, Geom::Vector<3,T> const & q1,
                              T & ratiop, T & ratioq)
  {
    Geom::Vector<3,T> pv(p1 - p0);
    Geom::Vector<3,T> qv(q1 - q0);
    Geom::Vector<3,T> n(unit(pv).cross(unit(qv)));
    // For very small angles, set normal to 0 to force the parallel
    // intersector code.
    T det = (n.length() < 0.01) ? 0.0 : determinant(pv, qv, n);
    if (det != 0.0)
    {
      Geom::Vector<3,T> pq(q0 - p0);
      ratiop =  determinant(pq, qv, n) / det;
      ratioq = -determinant(pv, pq, n) / det;
      return;
    }

    lineLineParallelPoint(p0, p1, q0, q1, pv, qv, ratiop, ratioq);
  }

  // Instantiate
  template void lineLineIntersectPnt_t(
                      Vector3d const & p0, Vector3d const & p1,
                      Vector3d const & q0, Vector3d const & q1,
                      double & ratiop, double & ratioq);


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
                            T tol, Cross_t<N,T> * cross)
  {
    double r;
    VectorNT pv(p1  - p0);
    VectorNT qv(pnt - p0);
    double pv_dot_pv = pv.dot(pv);
    double qv_dot_qv = qv.dot(qv);

    if (pv_dot_pv < tol*tol)
    {
      if (cross)
        *cross = Cross_t<N,T>(0.0);
      return qv_dot_qv;
    }

    if (cross)
      *cross = cross_t(pv, qv);

    T ratio = pv.dot(qv) / pv_dot_pv;

    if (ratio > 0.0 && ratio < 1.0)
    {
      r = qv_dot_qv - ratio*ratio * pv_dot_pv;
      return (r < 0.0) ? 0.0 : r;
    }
    else
    {
      if (ratio <= 0.0)
        return qv_dot_qv;
      else
      {
        VectorNT wk(qv - pv);
        return wk.dot(wk);
      }
    }
  }

  template double edgePointDistance2Pnt_t(Vector2d const & pnt,
                                          Vector2d const & p0,
                                          Vector2d const & p1,
                                          double tol,
                                          Cross_t<2,double> * cross);

  template double edgePointDistance2Pnt_t(Vector3d const & pnt,
                                          Vector3d const & p0,
                                          Vector3d const & p1,
                                          double tol,
                                          Cross_t<3,double> * cross);


  template < int N, typename T >
  void edgeEdgeIntersectPnt_t(VectorNT const &p0, VectorNT const &p1,
                              VectorNT const &q0, VectorNT const &q1,
                              T & ratiop, T & ratioq)
  {
    lineLineIntersectPnt_t(p0, p1, q0, q1, ratiop, ratioq);

    if (ratiop < 0.0) ratiop = 0.0;
    if (ratiop > 1.0) ratiop = 1.0;
    if (ratioq < 0.0) ratioq = 0.0;
    if (ratioq > 1.0) ratioq = 1.0;
  }

  // Instantiate
  template void edgeEdgeIntersectPnt_t(
                        Vector3d const &p0, Vector3d const &p1,
                        Vector3d const &q0, Vector3d const &q1,
                        double & ratiop, double & ratioq);


  template < int N, typename T >
  T edgeEdgeDistancePnt_t(VectorNT const & p0, VectorNT const & p1,
                          VectorNT const & q0, VectorNT const & q1,
                          T tol)
  {
#if 0
    std::cout << "v,1," << p0(0) << "," << p0(1) << "," << p0(2) << std::endl;
    std::cout << "v,2," << p1(0) << "," << p1(1) << "," << p1(2) << std::endl;
    std::cout << "v,3," << q0(0) << "," << q0(1) << "," << q0(2) << std::endl;
    std::cout << "v,4," << q1(0) << "," << q1(1) << "," << q1(2) << std::endl;

    VectorNT c(VectorNT(p1-p0).cross(VectorNT(q1-q0)));
    std::cout << "view," << c(0) << "," << c(1) << "," << c(2) << std::endl;
    std::cout << "repl" << std::endl;
#endif

    Cross_t<N,T> crs0, crs1;
    T dist2, minDist2 = std::numeric_limits < T >::max();

    dist2 = edgePointDistance2Pnt_t(q0, p0, p1, tol, &crs0);
    minDist2 = std::min(dist2, minDist2);
    dist2 = edgePointDistance2Pnt_t(q1, p0, p1, tol, &crs1);
    minDist2 = std::min(dist2, minDist2);
    int qWind = (dot_t(crs0, crs1) < 0.0) ? -1 : 1;

    dist2 = edgePointDistance2Pnt_t(p0, q0, q1, tol, &crs0);
    minDist2 = std::min(dist2, minDist2);
    dist2 = edgePointDistance2Pnt_t(p1, q0, q1, tol, &crs1);
    minDist2 = std::min(dist2, minDist2);
    int pWind = (dot_t(crs0, crs1) < 0.0) ? -1 : 1;

    if (qWind != -1 || pWind != -1)
      return std::sqrt(minDist2);

    return lineLineDistance(p0, p1, q0, q1);
  }

  template double edgeEdgeDistancePnt_t(
                     Vector3d const & p0, Vector3d const & p1,
                     Vector3d const & q0, Vector3d const & q1,
                     double tol);



#if 0
  // test for   bool edgeEdgeIntersect(Vector2d p0,
  //                                   Vector2d p1,
  //                                   Vector2d q0,
  //                                   Vector2d q1,
  //                                   double linearTol)

  // g++ -g -I ~/dev/cdna/star/base/src -I ~/dev/cdna/star/meshkernel/src test.cpp -L ~/dev/cdna/star/lib/linux-x86_64-2.6.1/gnu4.2-g/lib -l MKCore -l MKCommon -l MKLinalg -l MKMemory -l MKStorage -l Cad -o test

  // setenv LD_LIBRARY_PATH ~/dev/cdna/star/lib/linux-x86_64-2.6.1/gnu4.2-g/lib

#include <iostream>

// Spaces are needed here so the conversion to star-cd4's MK doesn't
// transform them.
# include "common/Vector.h"
# include "mk/core/EdgeCalculations2d.h"

  int main(int argc, char *argv[])
  {
    bool result;
    result = MKCore::edgeEdgeIntersect(Vector2d(0,0),
                                       Vector2d(1,1),
                                       Vector2d(0,1),
                                       Vector2d(1,0),
                                       0.0);
    std::cout << result << std::endl;

    result = MKCore::edgeEdgeIntersect(Vector2d(1,2),
                                       Vector2d(2,3),
                                       Vector2d(1,3),
                                       Vector2d(2,2),
                                       0.0);
    std::cout << result << std::endl;

    result = MKCore::edgeEdgeIntersect(Vector2d(0,0),
                                       Vector2d(1,1),
                                       Vector2d(0.5,0),
                                       Vector2d(0.5,1),
                                       0.0);
    std::cout << result << std::endl;

    result = MKCore::edgeEdgeIntersect(Vector2d(1,2),
                                       Vector2d(2,3),
                                       Vector2d(1.5,2),
                                       Vector2d(1.5,3),
                                       0.0);
    std::cout << result << std::endl;

    result = MKCore::edgeEdgeIntersect(Vector2d(0,0),
                                       Vector2d(1,1),
                                       Vector2d(1.5,0),
                                       Vector2d(1.5,3),
                                       0.0);
    std::cout << result << std::endl;

    result = MKCore::edgeEdgeIntersect(Vector2d(1,2),
                                       Vector2d(2,3),
                                       Vector2d(2.5,2),
                                       Vector2d(2.5,5),
                                       0.0);
    std::cout << result << std::endl;

    result = MKCore::edgeEdgeIntersect(Vector2d(0,0),
                                       Vector2d(1,1),
                                       Vector2d(1.5,1.4),
                                       Vector2d(1.5,1.6),
                                       0.0);
    std::cout << result << std::endl;

    result = MKCore::edgeEdgeIntersect(Vector2d(1,2),
                                       Vector2d(2,3),
                                       Vector2d(2.5,3.4),
                                       Vector2d(2.5,3.6),
                                       0.0);
    std::cout << result << std::endl;

    result = MKCore::edgeEdgeIntersect(Vector2d(0,0),
                                       Vector2d(1,1),
                                       Vector2d(0.5,0),
                                       Vector2d(0.5,0.5));
    std::cout << result << std::endl;

    result = MKCore::edgeEdgeIntersect(Vector2d(1,2),
                                       Vector2d(2,3),
                                       Vector2d(1.5,2),
                                       Vector2d(1.5,2.5));
    std::cout << result << std::endl;

    result = MKCore::edgeEdgeIntersect(Vector2d(0,0),
                                       Vector2d(1,1),
                                       Vector2d(0,-0.1),
                                       Vector2d(0,0.1));
    std::cout << result << std::endl;

    result = MKCore::edgeEdgeIntersect(Vector2d(1,2),
                                       Vector2d(2,3),
                                       Vector2d(1,-1.9),
                                       Vector2d(1,2.1));
    std::cout << result << std::endl;

    result = MKCore::edgeEdgeIntersect(Vector2d(0,0),
                                       Vector2d(1,1),
                                       Vector2d(1,0.9),
                                       Vector2d(1,1.1));
    std::cout << result << std::endl;

    result = MKCore::edgeEdgeIntersect(Vector2d(1,2),
                                       Vector2d(2,3),
                                       Vector2d(2,2.9),
                                       Vector2d(2,3.1));
    std::cout << result << std::endl;

    result = MKCore::edgeEdgeIntersect(
               Vector2d( 0.0038461574219799681, -0.051739125412470992),
               Vector2d(-0.024848446680010721,  -0.061454512211461948),
               Vector2d(-0.024779252138034447,  -0.053836001896839181),
               Vector2d(-0.024581501077207954,   0.035746681151061124),
               7.137672863998765e-05);
    std::cout << result << std::endl;
  }

  // The results should be:
  // 1
  // 1
  // 1
  // 1
  // 0
  // 0
  // 0
  // 0
  // 1
  // 1
  // 1
  // 1
  // 1
  // 1
  // 0
#endif


}
