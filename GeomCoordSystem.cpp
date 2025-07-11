#include "GeomCoordSystem.h"

#include <limits>

#include "GeomUtility.h"

namespace Geom
{
  // Find a complete 3D set of vectors that are orthonormal and righthanded.
  // The z axis is given in this->rot[2][0..2].  x and y are calculated.
  void CoordSystem::calcOrthoNormalSet()
  {
    rot[0] = Vector3d(1.0, 0.0, 0.0);
    rot[1] = Vector3d(0.0, 1.0, 0.0);
    rot[2].normalize();

    double dotx = rot[2](0);
    double doty = rot[2](1);

    // Use global x axis
    if (abs(dotx) < abs(doty))
    {
      rot[1] = rot[2].cross(rot[0]);
      rot[1].normalize();
      rot[0] = rot[1].cross(rot[2]);
      rot[0].normalize();
    }

    // Use global y axis
    {
      rot[0] = rot[1].cross(rot[2]);
      rot[0].normalize();
      rot[1] = rot[2].cross(rot[0]);
      rot[1].normalize();
    }

#if 0
    // I'm using the scalar triple product to find whither the x or y axis
    // forms the greatest angle with the z axis in absolute terms and start
    // with the local rotations in term of that axis. The sign of the STP
    // is used to keep the system righthanded.
    double STPx = rot[1].cross(rot[2]).dot(rot[0]);
    double STPy = rot[2].cross(rot[0]).dot(rot[1]);

    // Use global x axis
    if (abs(STPx) > abs(STPy))
    {
      rot[0](0) = (STPx < 0.0) ? -1.0 : 1.0;
      rot[1] = rot[2].cross(rot[0]);
      rot[1].normalize();
      rot[0] = rot[1].cross(rot[2]);
      rot[0].normalize();
    }

    // Use global y axis
    {
      rot[1](1) = (STPy < 0.0) ? -1.0 : 1.0;
      rot[0] = rot[1].cross(rot[2]);
      rot[0].normalize();
      rot[1] = rot[2].cross(rot[0]);
      rot[1].normalize();
    }
#endif
  }


  void CSmultiply(double const (&a)[4][4],
                  double const (&b)[4][4], double (&c)[4][4])
      {
        int i, j, k;
        for (i=0; i<4; ++i)
          for (j=0; j<4; ++j)
          {
            c[i][j] = 0.0;
            for (k=0; k<4; ++k)
              c[i][j] += a[i][k] * b[k][j];
          }
      }


  // If L1coord = [A] (rhs) Gcoord and Lcoord = [B] (*this) * L1coord,
  // that is, [A] is the first transformation and [B] is the second
  // transformation, we want to calculate [C] from [A] and [B] so that
  // Lcoord = [C] Gcoord.
  CoordSystem CoordSystem::chain(CoordSystem const & rhs) const
  {
    int i, j;
    double a[4][4], b[4][4], c[4][4];

    // untranspose matrices.
    for (i=0; i<3; ++i)
    {
      for (j=0; j<3; ++j)
      {
        a[i][j] = rhs.rot[j](i);
        b[i][j] = this->rot[j](i);
      }
      a[i][3] = rhs.origin(i);
      b[i][3] = this->origin(i);
      b[3][i] = a[3][i] = 0.0;
    }
    a[3][3] = b[3][3] = 1.0;


    // multiply.
    CSmultiply(a, b, c);


    // transpose matrix.
    Vector3d _origin, xaxis, yaxis, zaxis;
    for (i=0; i<3; ++i)
    {
      xaxis(i) = c[i][0];
      yaxis(i) = c[i][1];
      zaxis(i) = c[i][2];
      _origin(i) = c[i][3];
    }

    return CoordSystem(_origin, xaxis, yaxis, zaxis);
  }
}
