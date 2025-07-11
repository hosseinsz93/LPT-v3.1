#include "GeomXFaceCalculations.h"

#include "GeomBoundingBox.h"
#include "GeomXTol.h"
#include "GeomXUtilities.h"

#include "GeomXEdgeCalculations.h"
// #include "GeomXGenericVertices_t.h"


namespace GeomX
{
#ifndef NDEBUG
// #define SIDED_DEBUG
#endif

  // Is the point inside or outside the given face (polygon).  This routine
  // works in 2D since this operation only makes sense in that system.
  // Dimensions higher than 2 are ignored.

  // C is a container class that contains VectorNT.  It has an operator[]
  // to access elements of the array and a size() procedure that returns
  // how many elements are in the array.  N must be 2 or greater.  T
  // must be float or double.
  template < typename C, int N, typename T >
  Sidedness pointSideofFace(C const & argPoly,
                            VectorNT const & test, double linearTol)
  {
#ifndef TOL
#  define TOL 0.00001
#  define MACHINE_TOL 1.0e-30
#  define UNDEFTOL
#endif


#ifdef SIDED_DEBUG
    int testMin = 0, testMax = 1000000;
    static int cnt = 0;
    {
      ++cnt;
      std::cout << "\npointSideofFace " << cnt << std::endl;

      if (cnt >= testMin && cnt <= testMax)
      {
        unsigned i, j;
        std::cout << "dbase,clear" << std::endl;
        std::cout << "y" << std::endl;
        std::cout << "ctab,2,line,8" << std::endl;
        std::cout << "plty,norm" << std::endl;

        for (i=1; i<=argPoly.size(); ++i)
          std::cout << "v," << i << "," << argPoly[i-1](0)
                                 << "," << argPoly[i-1](1) << std::endl;
        std::cout << "v," << i << "," << test(0)
                               << "," << test(1) << std::endl;
        std::cout << "vset,news,vlist," << i << std::endl;

        for (i=1, j=2; i<=argPoly.size(); ++i, ++j)
        {
          if (j > argPoly.size()) j = 1;
          std::cout << "c," << i << "," << j <<std::endl;
        }
        std::cout << "cset,all" << std::endl;
        std::cout << "cdis,on,vert" << std::endl;
        std::cout << "cplo" << std::endl;
      }
    }
#endif


    VectorNT const * poly = &argPoly[0];
    VectorNT const * polyEnd = &argPoly[argPoly.size()];
    VectorNT xmin, xmax;

/*
 *  This algorithm works by counting the number of times a test line crosses
 *  the polygon edge.  The test line is a line from (x,y) to (-infinity,y).
 *
 *  Find the bounding box for the polygon.  If the test point is outside the
 *  box, it is outside the polygon.  If there is no linearTol, it will also
 *  be calculated here.
 */
    VectorNT const * pntc;
    {
      Geom::BoundingBox bb;
      for (pntc=poly; pntc<polyEnd; ++pntc)
        bb.expand(zeroZ((*pntc)));

      if (linearTol <= 0.0)
        linearTol = std::max(TOL * bb.longestSideLength(), MACHINE_TOL);

      if (Tol(test(0), linearTol) < bb.min()(0) ||
          Tol(test(0), linearTol) > bb.max()(0) ||
          Tol(test(1), linearTol) < bb.min()(1) ||
          Tol(test(1), linearTol) > bb.max()(1))
      {
#ifdef SIDED_DEBUG
        if (cnt >= testMin && cnt <= testMax)
          std::cout << "OUTSIDE\n" << std::endl;
#endif
        return OUTSIDE;
      }
    }


    // Find a point which does not lay on the test line.  This is so that
    // the code to check crossings when the polygon segment is on the test
    // line has a point off the line in order to see if the test line
    // crosses from going through the boundary.
    for (pntc=poly; pntc<polyEnd; ++pntc)
    {
      if (Tol((*pntc)(1), linearTol) != test(1)) break;
      if (Tol((*pntc)(0), linearTol) == test(0))
      {
#ifdef SIDED_DEBUG
        if (cnt >= testMin && cnt <= testMax)
          std::cout << "BOUNDARY\n" << std::endl;
#endif
        return BOUNDARY;
      }
      if (Tol((*pntc)(0), linearTol) > test(0)) break;
      // If you get to here, you know the polygon point is on the test line
      // (remember, it goes from (x,y) to (-infinity,y)).
    }
    if (pntc >= polyEnd)
    {
      // If you get to here, all the polygon points lay on the test line
      // and the test point must not be in the polygon.
#ifdef SIDED_DEBUG
      if (cnt >= testMin && cnt <= testMax)
        std::cout << "OUTSIDE\n" << std::endl;
#endif
      return OUTSIDE;
    }


    // Visit all the polygon edges starting with the point found above.
    unsigned i, size=argPoly.size(), cross = 0;
    VectorNT const * pntn;
    VectorNT const * beforeBound = 0;
    double xminBound = 0.0, xmaxBound = 0.0;
    for (i=0, /*pntc=pntc,*/ pntn=pntc, ++pntn; i<size; ++i, pntc=pntn, ++pntn)
    {
      if (pntn >= polyEnd) pntn = poly;

      // If the current polygon segment end point is on the test line,
      // save information so you can figure out the crossing when you
      // come to an end point that isn't on the test line.  Keep track
      // of the span of x to see if the point is on the boundary.
      if (Tol((*pntn)(1), linearTol) == test(1))
      {
        if (beforeBound)
        {
          xminBound = std::min(xminBound, (*pntn)(0));
          xmaxBound = std::max(xmaxBound, (*pntn)(0));
          if (Tol(test(0), linearTol) >= xminBound &&
              Tol(test(0), linearTol) <= xmaxBound)
          {
#ifdef SIDED_DEBUG
            if (cnt >= testMin && cnt <= testMax)
              std::cout << "BOUNDARY\n" << std::endl;
#endif
            return BOUNDARY;
          }
          continue;
        }

        if (Tol((*pntn)(0), linearTol) == test(0))
        {
#ifdef SIDED_DEBUG
          if (cnt >= testMin && cnt <= testMax)
            std::cout << "BOUNDARY\n" << std::endl;
#endif
          return BOUNDARY;
        }

        beforeBound = pntc;
        xminBound = xmaxBound = (*pntn)(0);
        continue;
      }


      // If we've just come off the test line, check to see if transvering
      // the boundary results in a crossing;
      if (beforeBound)
      {
        // boundary is on "wrong" side, there cannot be a crossing.  So only
        // do check if it is on the correct side.
        if (Tol(test(0), linearTol) > xmaxBound)
        {
          double sign = ((*beforeBound)(1) - test(1)) * ((*pntn)(1) - test(1));
          if (sign < 0.0) ++cross;
        }

        beforeBound = 0;
        continue;
      }


      // Now we take a look at the polygon segment and see if there is a
      // crossing.  First do the gross tests (bounding box) then the more
      // "refined" tests if necessary.  The order of the test is important
      // here.  Don't change them.
      Geom::BoundingBox bb(*pntc);
      bb.expand(*pntn);
      if (Tol(test(0), linearTol) < bb.min()(0)) continue;
      if (Tol(test(1), linearTol) < bb.min()(1)) continue;
      if (Tol(test(1), linearTol) > bb.max()(1)) continue;
      if (Tol(test(0), linearTol) > bb.max()(0))
      {
        ++cross;
        continue;
      }

      // The refined test.  Calculate the location of the intersection of the
      // test line and the current line segment.
      double invSlope = ((*pntn)(0) - (*pntc)(0)) / ((*pntn)(1) - (*pntc)(1));
      double x = (*pntc)(0) - invSlope * ((*pntc)(1) - test(1));
      if (Tol(test(0), linearTol) == x)
      {
#ifdef SIDED_DEBUG
        if (cnt >= testMin && cnt <= testMax)
          std::cout << "BOUNDARY\n" << std::endl;
#endif
        return BOUNDARY;
      }
      if (Tol(test(0), linearTol) > x)
        ++cross;
    }

    if (cross % 2 != 0)
    {
#ifdef SIDED_DEBUG
      if (cnt >= testMin && cnt <= testMax)
        std::cout << "INSIDE\n" << std::endl;
#endif
      return INSIDE;
    }

#ifdef SIDED_DEBUG
    if (cnt >= testMin && cnt <= testMax)
      std::cout << "OUTSIDE\n" << std::endl;
#endif
    return OUTSIDE;

#ifdef UNDEFTOL
#  undef TOL
#  undef MACHINE_TOL
#  undef UNDEFTOL
#endif
  }



  // C is a container class that contains VectorNT.  It has an operator[]
  // to access elements of the array and a size() procedure that returns
  // how many elements are in the array.  N must be 3 and T must be double.
  template < typename C >
  bool pointInsideOnFace(C const & coords, Vector3d p)
  {
    Sidedness sidedness = pointSideofFace(coords, p);
    return sidedness != OUTSIDE;
  }


  template < int N, typename T >
  double pointFaceDistance2(VectorNT p,
                            std::vector < VectorNT > const & fcoords,
                            double linearTol)
  {
    if (pointInsideOnFace(fcoords, p))
      return 0.0;

    typename std::vector < VectorNT >::const_iterator i, j;
    double minDist2 = std::numeric_limits<double>::max();
    for (i=fcoords.begin(), j=i+1; i!=fcoords.end(); ++i, ++j)
    {
      if (j == fcoords.end()) j=fcoords.begin();
      minDist2 = std::min(minDist2,
                          edgePointDistance2Pnt_t(p, *i, *j, linearTol));
    }
    return minDist2;
  }
}
