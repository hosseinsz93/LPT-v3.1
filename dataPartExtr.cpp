static char help[] = "Testing programming!";

#include <algorithm>    // std::sort
#include <array>
#include <fstream>
#include <limits>
#include <math.h>
#include <memory>
#include <regex>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "petscda.h"
#include "petscts.h"
#include "petscpc.h"
#include "petscsnes.h"

#include "DragCalcDebug.h"
#include "GeomRtree.h"
#include "GeomSubset.h"
#include "GeomXTol.h"
#include "Id.h"
#include "SpatialIndexGrid.h"
#include "timer.h"

  
#define NEWMETRIC

#include "TECIO.h"

// Needed to get tracking.o to link correctly
PetscReal ptc_dt;
PetscInt ptc_ftn;
PetscInt ptc_int;
PetscReal cl;


PetscInt ti, block_number, Flux_in;
int binary_input=0;
int xyz_input=0;
PetscInt tis, tie, tsteps=5;
PetscReal angle;
int nv_once=0;
int onlyV=0;
int k_average=0;
int j_average=0;
int i_average=0;
int ik_average=0;
int ikc_average=0;      // conditional spatial averaging in ik directions (channel flow)
int reynolds=0; // 1: contravariant reynolds stress

int i_begin, i_end;
int j_begin, j_end;
int k_begin, k_end;

int pcr=0;
int avg=0, rans=0, rans_output=0, levelset=0, conv_diff=0, sandwave=0;
int vc = 1;
int width_ave_instant = 0;

int cs=0;
int i_periodic=0;
int j_periodic=0;
int k_periodic=0;
int kk_periodic=0;
int averaging_option=0;
int pi=-1, pk=-1;
int shear=0;

char prefix[256];

double gravity_x=0, gravity_y=0, gravity_z=0;

//int l, m, n;
/* Symmetric matrix A -> eigenvectors in columns of V, corresponding eigenvalues in d. */
void eigen_decomposition(double A[3][3], double V[3][3], double d[3]);


typedef struct {
    PassiveScalar u, v, w, p;
} PassiveField;

typedef struct {
    PetscScalar u, v, w;
} Field;

#include "cmpnts.h"

// typedef struct {
//     PetscScalar x, y, z;
// } Cmpnts;

typedef struct {
    PetscScalar x, y;
} Cmpnts2;

typedef struct {
    PassiveScalar csi[3], eta[3], zet[3], aj;
} Metrics;


typedef struct {
    PetscInt      nbnumber;
    PetscInt      n_v, n_elmt;    // number of vertices and number of elements
    PetscInt      *nv1, *nv2, *nv3;       // Node index of each triangle
    PetscReal     *nf_x, *nf_y, *nf_z;    // Normal direction
    PetscReal     *x_bp, *y_bp, *z_bp;    // Coordinates of IBM surface nodes
    PetscReal     *x_bp0, *y_bp0, *z_bp0;
    PetscReal     *x_bp_o, *y_bp_o, *z_bp_o;
/*   PetscReal  x_bp_in[101][3270], y_bp_in[101][3270], z_bp_in[101][3270]; */
    Cmpnts        *u, *uold;

} IBMNodes;

#include "UserCtx.h"

// Predefinitions of subroutines defined later.
PetscErrorCode ReadCoordinates(std::shared_ptr< UserCtx > & user);
PetscErrorCode QCriteria(std::shared_ptr< UserCtx > & user);
PetscErrorCode Velocity_Magnitude(std::shared_ptr< UserCtx > & user);
PetscErrorCode Velocity_Magnitude_4(std::shared_ptr<UserCtx> & user, PetscInt);
PetscErrorCode Lambda2(std::shared_ptr< UserCtx > & user);
PetscErrorCode FormMetrics(std::shared_ptr< UserCtx > & user);
void Calc_avg_shear_stress(std::shared_ptr< UserCtx > & user);
void outputUandFTLE(std::shared_ptr< UserCtx > & user, int ti);

namespace LOCAL
{
  void genRtree(std::shared_ptr< UserCtx > & user);
}
void calcPartCount(std::shared_ptr< UserCtx > & user, int ti);
void calcXCrossing(std::shared_ptr< UserCtx > & user, int ti);
void crossXPartOut(std::shared_ptr< UserCtx > & user,
                   std::vector<int> & sortedpId);
void crossXContOut(std::shared_ptr< UserCtx > & user,
                   std::vector<int> & sortedpId);
void clearCount(std::shared_ptr< UserCtx > & user);
void calcCount(std::shared_ptr< UserCtx > & user, int ti);
void concOutput(std::shared_ptr< UserCtx > & user, int ti, int fileCnt);
void calcZMaxCount(std::shared_ptr< UserCtx > & user, int ti);


template < typename T > T sqr(T const & x) { return x*x; }

// Interpolate data that is not exactly on the mesh onto the mesh.
// I'm using 1/d, where d is the distance of the known value to the
// mesh point, as the distribution weight.

// axis identifies which i,j,k indices are on the plane([0] and [1]) and
// which index is off the plane [2]. The calculations are done on the defined
// plane. Use 0 for i, 1 for j and 2 for k.

// i, j, k are the indices closeset to x,y,z
// x, y, z are the coordinates of the value location
// value is the value at x,y,z

class InterpolateVectorData
{
  public:
    InterpolateVectorData(int xs, int xe, int ys, int ye, int zs, int ze);


    void operator()(int const (&axis)[3],
                    int const i, int const j, int const k,
                    Vector3d const & c, Cmpnts ***coor,
                    double val, PetscReal ***varr);
    
  private:
    int se[3][2];
};


// WRO - class to write paraview readable vtk files.

int paraview = 0;
// Is the paraview file going to be written in ascii or binary?
bool paraviewBinFile = false;

class Paraview
{
  public:
  Paraview(Paraview const &);
  Paraview & operator=(Paraview const &);

  ~Paraview()
  {
    if (_ptr && f)
      fclose(f);
  }

  static void init(char const * filename)
  {
    if (_ptr)
    {
      printf("Paraview class is already intialized.\n");
      exit(1);
    }

    _ptr = new Paraview(paraviewBinFile);

    _ptr->f = fopen(filename, "w");
    if (!_ptr->f)
    {
      printf("Couldn't open file %s.\n", filename);
      exit(1);
    }

    fprintf(_ptr->f, "# vtk DataFile Version 3.0\n");
    fprintf(_ptr->f,
            "Describe the data and include any other pertinent information\n");
    fprintf(_ptr->f, _ptr->binary ? "BINARY\n" : "ASCII\n");
  }

  static Paraview * ptr()
  {
    if (_ptr)
      return _ptr;

    printf("Class Paraview has not been initialized yet!");
    exit(1);
  }

  void outputGrid(PetscInt xs, PetscInt xe,
                  PetscInt ys, PetscInt ye,
                  PetscInt zs, PetscInt ze, Cmpnts ***coor)
  {
    printf("writing %s file.\n", binary ? "binary" : "ascii");
    fprintf(f, "DATASET STRUCTURED_GRID\n");
    PetscInt xn = xe - xs - 1;
    PetscInt yn = ye - ys - 1;
    PetscInt zn = ze - zs - 1;
    fprintf(f, "DIMENSIONS %d %d %d\n", xn, yn, zn);
    fprintf(f, "POINTS %d float\n", xn * yn * zn);

    for (PetscInt k=zs; k<ze-1; k++)
      for (PetscInt j=ys; j<ye-1; j++)
        for (PetscInt i=xs; i<xe-1; i++)
        {
          if (binary)
          {
            if (fwrite(&coor[k][j][i].x, sizeof(float), 3, f) != 3)
            {
              printf("Error writing binary grid to file.\n");
              exit(1);
            }
          }
          else
          {
            if (fprintf(f, "%.7f %.7f %.7f\n",
                        coor[k][j][i].x, coor[k][j][i].y, coor[k][j][i].z) < 0)
            {
              printf("Error writing ascii grid to file.\n");
              exit(1);
            }
          }
        }
  }

  void outputData(char const * desc, float * x, int cnt){
    fprintf(f, "POINT_DATA %d\n", cnt);
    fprintf(f, "SCALARS %s float 1\n", desc);
    fprintf(f, "LOOKUP_TABLE default\n");

    if (binary)
    {
      if (fwrite(x, sizeof(float), cnt, f) != cnt)
      {
        printf("Error writing binary %d to file.\n", desc);
        exit(1);
      }
    }
    else
    {
      for (int i=0; i<cnt; ++i)
      {
        if (fprintf(f, "%.7f\n", x[i]) < 0)
        {
          printf("Error writing ascii grid to file.\n");
          exit(1);
        }
      }
    }
  }


  private:
  static Paraview * _ptr;

  FILE * f;
  bool binary;


  Paraview(bool _binary = true)                            // constructor is private
    : f(NULL)
    , binary(_binary)
  { };
};

Paraview * Paraview::_ptr = NULL;



int file_exist(char *str)
{
  int r=0;

  /*if(!my_rank)*/ {
    FILE *fp=fopen(str, "r");
    if(!fp) {
      r=0;
      printf("\nFILE !!! %s does not exist !!!\n", str);
    }
    else {
      fclose(fp);
      r=1;
    }
  }
  MPI_Bcast(&r, 1, MPI_INT, 0, PETSC_COMM_WORLD);
  return r;
};

void Calculate_Covariant_metrics(double g[3][3], double G[3][3])
{
  /*
    | csi.x  csi.y csi.z |-1                | x.csi  x.eta x.zet |
    | eta.x eta.y eta.z |    =      | y.csi   y.eta  y.zet |
    | zet.x zet.y zet.z |           | z.csi  z.eta z.zet |

  */
  const double a11=g[0][0], a12=g[0][1], a13=g[0][2];
  const double a21=g[1][0], a22=g[1][1], a23=g[1][2];
  const double a31=g[2][0], a32=g[2][1], a33=g[2][2];

  double det= a11*(a33*a22-a32*a23)-a21*(a33*a12-a32*a13)+a31*(a23*a12-a22*a13);

  G[0][0] = (a33*a22-a32*a23)/det,        G[0][1] = - (a33*a12-a32*a13)/det,      G[0][2] = (a23*a12-a22*a13)/det;
  G[1][0] = -(a33*a21-a31*a23)/det, G[1][1] = (a33*a11-a31*a13)/det,      G[1][2] = - (a23*a11-a21*a13)/det;
  G[2][0] = (a32*a21-a31*a22)/det,        G[2][1] = - (a32*a11-a31*a12)/det,      G[2][2] = (a22*a11-a21*a12)/det;
};

void Calculate_normal(Cmpnts csi, Cmpnts eta, Cmpnts zet, double ni[3], double nj[3], double nk[3])
{
  double g[3][3];
  double G[3][3];

  g[0][0]=csi.x, g[0][1]=csi.y, g[0][2]=csi.z;
  g[1][0]=eta.x, g[1][1]=eta.y, g[1][2]=eta.z;
  g[2][0]=zet.x, g[2][1]=zet.y, g[2][2]=zet.z;

  Calculate_Covariant_metrics(g, G);
  double xcsi=G[0][0], ycsi=G[1][0], zcsi=G[2][0];
  double xeta=G[0][1], yeta=G[1][1], zeta=G[2][1];
  double xzet=G[0][2], yzet=G[1][2], zzet=G[2][2];

  double nx_i = xcsi, ny_i = ycsi, nz_i = zcsi;
  double nx_j = xeta, ny_j = yeta, nz_j = zeta;
  double nx_k = xzet, ny_k = yzet, nz_k = zzet;

  double sum_i=sqrt(nx_i*nx_i+ny_i*ny_i+nz_i*nz_i);
  double sum_j=sqrt(nx_j*nx_j+ny_j*ny_j+nz_j*nz_j);
  double sum_k=sqrt(nx_k*nx_k+ny_k*ny_k+nz_k*nz_k);

  nx_i /= sum_i, ny_i /= sum_i, nz_i /= sum_i;
  nx_j /= sum_j, ny_j /= sum_j, nz_j /= sum_j;
  nx_k /= sum_k, ny_k /= sum_k, nz_k /= sum_k;

  ni[0] = nx_i, ni[1] = ny_i, ni[2] = nz_i;
  nj[0] = nx_j, nj[1] = ny_j, nj[2] = nz_j;
  nk[0] = nx_k, nk[1] = ny_k, nk[2] = nz_k;
}

double Contravariant_Reynolds_stress(double uu, double uv, double uw, double vv, double vw, double ww,
                                     double csi0, double csi1, double csi2, double eta0, double eta1, double eta2)
{

  double A = uu*csi0*eta0 + vv*csi1*eta1 + ww*csi2*eta2 + uv * (csi0*eta1+csi1*eta0)      + uw * (csi0*eta2+csi2*eta0) + vw * (csi1*eta2+csi2*eta1);
  double B = sqrt(csi0*csi0+csi1*csi1+csi2*csi2)*sqrt(eta0*eta0+eta1*eta1+eta2*eta2);

  return A/B;
}

void Compute_du_center (int i, int j, int k,  int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert,
                        double *dudc, double *dvdc, double *dwdc,
                        double *dude, double *dvde, double *dwde,
                        double *dudz, double *dvdz, double *dwdz)
{
  if ((nvert[k][j][i+1])> 0.1 || (i==mx-2 && !i_periodic)) {
    *dudc = ( ucat[k][j][i].x - ucat[k][j][i-1].x );
    *dvdc = ( ucat[k][j][i].y - ucat[k][j][i-1].y );
    *dwdc = ( ucat[k][j][i].z - ucat[k][j][i-1].z );
  }
  else if ((nvert[k][j][i-1])> 0.1 || (i==1 && !i_periodic)) {
    *dudc = ( ucat[k][j][i+1].x - ucat[k][j][i].x );
    *dvdc = ( ucat[k][j][i+1].y - ucat[k][j][i].y );
    *dwdc = ( ucat[k][j][i+1].z - ucat[k][j][i].z );
  }
  else {
    if(i_periodic && i==1) {
      *dudc = ( ucat[k][j][i+1].x - ucat[k][j][mx-2].x ) * 0.5;
      *dvdc = ( ucat[k][j][i+1].y - ucat[k][j][mx-2].y ) * 0.5;
      *dwdc = ( ucat[k][j][i+1].z - ucat[k][j][mx-2].z ) * 0.5;
    }
    else if(i_periodic && i==mx-2) {
      *dudc = ( ucat[k][j][1].x - ucat[k][j][i-1].x ) * 0.5;
      *dvdc = ( ucat[k][j][1].y - ucat[k][j][i-1].y ) * 0.5;
      *dwdc = ( ucat[k][j][1].z - ucat[k][j][i-1].z ) * 0.5;
    }
    else {
      *dudc = ( ucat[k][j][i+1].x - ucat[k][j][i-1].x ) * 0.5;
      *dvdc = ( ucat[k][j][i+1].y - ucat[k][j][i-1].y ) * 0.5;
      *dwdc = ( ucat[k][j][i+1].z - ucat[k][j][i-1].z ) * 0.5;
    }
  }

  if ((nvert[k][j+1][i])> 0.1 || (j==my-2 && !j_periodic)) {
    *dude = ( ucat[k][j][i].x - ucat[k][j-1][i].x );
    *dvde = ( ucat[k][j][i].y - ucat[k][j-1][i].y );
    *dwde = ( ucat[k][j][i].z - ucat[k][j-1][i].z );
  }
  else if ((nvert[k][j-1][i])> 0.1 || (j==1 && !j_periodic)) {
    *dude = ( ucat[k][j+1][i].x - ucat[k][j][i].x );
    *dvde = ( ucat[k][j+1][i].y - ucat[k][j][i].y );
    *dwde = ( ucat[k][j+1][i].z - ucat[k][j][i].z );
  }
  else {
    if(j_periodic && j==1) {
      *dude = ( ucat[k][j+1][i].x - ucat[k][my-2][i].x ) * 0.5;
      *dvde = ( ucat[k][j+1][i].y - ucat[k][my-2][i].y ) * 0.5;
      *dwde = ( ucat[k][j+1][i].z - ucat[k][my-2][i].z ) * 0.5;
    }
    else if(j_periodic && j==my-2) {
      *dude = ( ucat[k][1][i].x - ucat[k][j-1][i].x ) * 0.5;
      *dvde = ( ucat[k][1][i].y - ucat[k][j-1][i].y ) * 0.5;
      *dwde = ( ucat[k][1][i].z - ucat[k][j-1][i].z ) * 0.5;
    }
    else {
      *dude = ( ucat[k][j+1][i].x - ucat[k][j-1][i].x ) * 0.5;
      *dvde = ( ucat[k][j+1][i].y - ucat[k][j-1][i].y ) * 0.5;
      *dwde = ( ucat[k][j+1][i].z - ucat[k][j-1][i].z ) * 0.5;
    }
  }

  if ((nvert[k+1][j][i])> 0.1 || (k==mz-2 && !k_periodic)) {
    *dudz = ( ucat[k][j][i].x - ucat[k-1][j][i].x );
    *dvdz = ( ucat[k][j][i].y - ucat[k-1][j][i].y );
    *dwdz = ( ucat[k][j][i].z - ucat[k-1][j][i].z );
  }
  else if ((nvert[k-1][j][i])> 0.1 || (k==1 && !k_periodic)) {
    *dudz = ( ucat[k+1][j][i].x - ucat[k][j][i].x );
    *dvdz = ( ucat[k+1][j][i].y - ucat[k][j][i].y );
    *dwdz = ( ucat[k+1][j][i].z - ucat[k][j][i].z );
  }
  else {
    if(k_periodic && k==1) {
      *dudz = ( ucat[k+1][j][i].x - ucat[mz-2][j][i].x ) * 0.5;
      *dvdz = ( ucat[k+1][j][i].y - ucat[mz-2][j][i].y ) * 0.5;
      *dwdz = ( ucat[k+1][j][i].z - ucat[mz-2][j][i].z ) * 0.5;
    }
    else if(k_periodic && k==mz-2) {
      *dudz = ( ucat[1][j][i].x - ucat[k-1][j][i].x ) * 0.5;
      *dvdz = ( ucat[1][j][i].y - ucat[k-1][j][i].y ) * 0.5;
      *dwdz = ( ucat[1][j][i].z - ucat[k-1][j][i].z ) * 0.5;
    }
    else {
      *dudz = ( ucat[k+1][j][i].x - ucat[k-1][j][i].x ) * 0.5;
      *dvdz = ( ucat[k+1][j][i].y - ucat[k-1][j][i].y ) * 0.5;
      *dwdz = ( ucat[k+1][j][i].z - ucat[k-1][j][i].z ) * 0.5;
    }
  }
}

void Compute_dscalar_center (int i, int j, int k,  int mx, int my, int mz, PetscReal ***p, PetscReal ***nvert, double *dpdc, double *dpde, double *dpdz)
{
  if ((nvert[k][j][i+1])> 0.1 || (i==mx-2 && !i_periodic)) {
    *dpdc = ( p[k][j][i] - p[k][j][i-1] );
  }
  else if ((nvert[k][j][i-1])> 0.1 || (i==1 && !i_periodic)) {
    *dpdc = ( p[k][j][i+1] - p[k][j][i] );
  }
  else {
    if(i_periodic && i==1) {
      *dpdc = ( p[k][j][i+1] - p[k][j][mx-2] ) * 0.5;
    }
    else if(i_periodic && i==mx-2) {
      *dpdc = ( p[k][j][1] - p[k][j][i-1] ) * 0.5;
    }
    else {
      *dpdc = ( p[k][j][i+1] - p[k][j][i-1] ) * 0.5;
    }
  }

  if ((nvert[k][j+1][i])> 0.1 || (j==my-2 && !j_periodic)) {
    *dpde = ( p[k][j][i] - p[k][j-1][i] );
  }
  else if ((nvert[k][j-1][i])> 0.1 || (j==1 && !j_periodic)) {
    *dpde = ( p[k][j+1][i] - p[k][j][i] );
  }
  else {
    if(j_periodic && j==1) {
      *dpde = ( p[k][j+1][i] - p[k][my-2][i] ) * 0.5;
    }
    else if(j_periodic && j==my-2) {
      *dpde = ( p[k][1][i] - p[k][j-1][i] ) * 0.5;
    }
    else {
      *dpde = ( p[k][j+1][i] - p[k][j-1][i] ) * 0.5;
    }
  }

  if ((nvert[k+1][j][i])> 0.1 || (k==mz-2 && !k_periodic)) {
    *dpdz = ( p[k][j][i] - p[k-1][j][i] );
  }
  else if ((nvert[k-1][j][i])> 0.1 || (k==1 && !k_periodic)) {
    *dpdz = ( p[k+1][j][i] - p[k][j][i] );
  }
  else {
    if(k_periodic && k==1) {
      *dpdz = ( p[k+1][j][i] - p[mz-2][j][i] ) * 0.5;
    }
    else if(k_periodic && k==mz-2) {
      *dpdz = ( p[1][j][i] - p[k-1][j][i] ) * 0.5;
    }
    else {
      *dpdz = ( p[k+1][j][i] - p[k-1][j][i] ) * 0.5;
    }
  }
}

void Compute_du_dxyz (  double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc,
                        double dudc, double dvdc, double dwdc, double dude, double dvde, double dwde, double dudz, double dvdz, double dwdz,
                        double *du_dx, double *dv_dx, double *dw_dx, double *du_dy, double *dv_dy, double *dw_dy, double *du_dz, double *dv_dz, double *dw_dz )
{
  *du_dx = (dudc * csi0 + dude * eta0 + dudz * zet0) * ajc;
  *du_dy = (dudc * csi1 + dude * eta1 + dudz * zet1) * ajc;
  *du_dz = (dudc * csi2 + dude * eta2 + dudz * zet2) * ajc;
  *dv_dx = (dvdc * csi0 + dvde * eta0 + dvdz * zet0) * ajc;
  *dv_dy = (dvdc * csi1 + dvde * eta1 + dvdz * zet1) * ajc;
  *dv_dz = (dvdc * csi2 + dvde * eta2 + dvdz * zet2) * ajc;
  *dw_dx = (dwdc * csi0 + dwde * eta0 + dwdz * zet0) * ajc;
  *dw_dy = (dwdc * csi1 + dwde * eta1 + dwdz * zet1) * ajc;
  *dw_dz = (dwdc * csi2 + dwde * eta2 + dwdz * zet2) * ajc;
};

void Compute_dscalar_dxyz (double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc,
                           double dpdc, double dpde, double dpdz, double *dp_dx, double *dp_dy, double *dp_dz)
{
  *dp_dx = (dpdc * csi0 + dpde * eta0 + dpdz * zet0) * ajc;
  *dp_dy = (dpdc * csi1 + dpde * eta1 + dpdz * zet1) * ajc;
  *dp_dz = (dpdc * csi2 + dpde * eta2 + dpdz * zet2) * ajc;
};


void IKavg(float *x, int xs, int xe, int ys, int ye, int zs, int ze, int mx, int my, int mz)
{
  int i, j, k;

  for (j=ys; j<ye-2; j++) {
    double iksum=0;
    int count=0;
    for (i=xs; i<xe-2; i++)
      for (k=zs; k<ze-2; k++) {
        iksum += x[k * (mx-2)*(my-2) + j*(mx-2) + i];
        count++;
      }
    double ikavg = iksum/(double)count;
    for (i=xs; i<xe-2; i++)
      for (k=zs; k<ze-2; k++)         x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ikavg;
  }
}

/*
  pi, pk : # of grid points correcsponding to the period
  conditional averaging
*/
void IKavg_c (float *x, int xs, int xe, int ys, int ye, int zs, int ze, int mx, int my, int mz)
{
  //int i, j, k;

  if(pi<=0) pi = (xe-xs-2); // no averaging
  if(pk<=0) pk = (ze-zs-2); // no averaging

  int ni = (xe-xs-2) / pi;
  int nk = (ze-zs-2) / pk;

  std::vector< std::vector<float> > iksum (pk);

  for(int k=0; k<pk; k++) iksum[k].resize(pi);

  for (int j=ys; j<ye-2; j++) {

    for(int k=0; k<pk; k++) std::fill( iksum[k].begin(), iksum[k].end(), 0.0 );

    for (int i=xs; i<xe-2; i++)
      for (int k=zs; k<ze-2; k++) {
        iksum [ ( k - zs ) % pk ] [ ( i - xs ) % pi] += x [k * (mx-2)*(my-2) + j*(mx-2) + i];
      }

    for (int i=xs; i<xe-2; i++)
      for (int k=zs; k<ze-2; k++) x [k * (mx-2)*(my-2) + j*(mx-2) + i] = iksum [ ( k - zs ) % pk ] [ ( i - xs ) % pi ] / (ni*nk);
  }
}


void Kavg(float *x, int xs, int xe, int ys, int ye, int zs, int ze, int mx, int my, int mz)
{
  int i, j, k;

  for (j=ys; j<ye-2; j++)
    for (i=xs; i<xe-2; i++) {
      double ksum=0;
      int count=0;
      for (k=zs; k<ze-2; k++) {
        ksum += x[k * (mx-2)*(my-2) + j*(mx-2) + i];
        count++;
      }
      double kavg = ksum/(double)count;
      for (k=zs; k<ze-2; k++) {
        x[k * (mx-2)*(my-2) + j*(mx-2) + i] = kavg;
      }
    }

}

void Javg(float *x, int xs, int xe, int ys, int ye, int zs, int ze, int mx, int my, int mz)
{
  int i, j, k;

  for (k=zs; k<ze-2; k++)
    for (i=xs; i<xe-2; i++) {
      double jsum=0;
      int count=0;
      for (j=ys; j<ye-2; j++) {
        jsum += x[k * (mx-2)*(my-2) + j*(mx-2) + i];
        count++;
      }
      double javg = jsum/(double)count;
      for (j=ys; j<ye-2; j++) {
        x[k * (mx-2)*(my-2) + j*(mx-2) + i] = javg;
      }
    }

}

void Iavg(float *x, int xs, int xe, int ys, int ye, int zs, int ze, int mx, int my, int mz)
{
  int i, j, k;

  for (k=zs; k<ze-2; k++) {
    for (j=ys; j<ye-2; j++) {
      double isum=0;
      int count=0;
      for (i=xs; i<xe-2; i++) {
        isum += x[k * (mx-2)*(my-2) + j*(mx-2) + i];
        count++;
      }
      double iavg = isum/(double)count;
      for (i=xs; i<xe-2; i++) {
        x[k * (mx-2)*(my-2) + j*(mx-2) + i] = iavg;
      }
    }
  }
}

void Iavg(Cmpnts ***x, int xs, int xe, int ys, int ye, int zs, int ze, int mx, int my, int mz)
{
  int i, j, k;

  for (k=zs; k<ze-2; k++)
    for (j=ys; j<ye-2; j++) {
      Cmpnts isum, iavg;
      isum.x = isum.y = isum.z = 0;

      int count=0;
      for (i=xs; i<xe-2; i++) {
        isum.x += x[k+1][j+1][i+1].x;
        isum.y += x[k+1][j+1][i+1].y;
        isum.z += x[k+1][j+1][i+1].z;
        count++;
      }

      iavg.x = isum.x /(double)count;
      iavg.y = isum.y /(double)count;
      iavg.z = isum.z /(double)count;

      for (i=xs; i<xe-2; i++) {
        x[k+1][j+1][i+1] = iavg;
      }
    }
}

void Iavg(PetscReal ***x, int xs, int xe, int ys, int ye, int zs, int ze, int mx, int my, int mz)
{
  int i, j, k;

  for (k=zs; k<ze-2; k++)
    for (j=ys; j<ye-2; j++) {
      double isum, iavg;
      isum = 0;

      int count=0;
      for (i=xs; i<xe-2; i++) {
        isum += x[k+1][j+1][i+1];
        count++;
      }
      iavg = isum /(double)count;

      for (i=xs; i<xe-2; i++) {
        x[k+1][j+1][i+1] = iavg;
      }
    }
}


// Find max value and output it.  WRO
class StoreMaxEntry
{
  public:
  StoreMaxEntry(int ia, int ja, int ka, double xa)
    : i(ia), j(ja), k(ka), x(xa)
  { }

  int i, j, k;
  double x;
};

template < int MaxSize = 10 > class StoreMax
{
  public:
  StoreMax(std::string const & vara) : var(vara) { }

  void checkEntry(int i, int j, int k, double x)
  {
    // This section of code finds the idx of the number in the table
    // that is >= than x.
    int idx;
    for (idx=0; idx<tab.size(); ++idx)
    {
      if (tab[idx].x < x)
        continue;

      break;
    }

    // x is not greater than anything in table, so don't add
    // unless there's nothing there anyway.
    if (idx == 0 && tab.size() != 0)
      return;

    // If table is not full, insert item at idx.
    if (tab.size() < MaxSize)
    {
      tab.insert(tab.begin()+idx, StoreMaxEntry(i, j, k, x));
      return;
    }


    // If table is full, insert new entry at idx
    // and remove item 0;
    for (int idx1=0; idx1<idx-1; ++idx1)
      tab[idx1] = tab[idx1+1];
    tab[idx-1] = StoreMaxEntry(i, j, k, x);
  }

  // Output the results.
  void output() const
  {
    printf("\nMaximum values for %s.\n", var.c_str());
    for (int idx=tab.size()-1; idx!=-1; --idx)
      printf("%f - %d, %d, %d\n", tab[idx].x, tab[idx].i, tab[idx].j, tab[idx].k);
    printf("\n");
  }

  std::string var;
  std::vector < StoreMaxEntry > tab;
};

class GridIndices
{
  public:
    PetscInt xs; PetscInt xe;
    PetscInt ys; PetscInt ye;
    PetscInt zs; PetscInt ze;
    PetscInt mx; PetscInt my;
    PetscInt mz;
    INTEGER4 IMax, JMax, KMax;

    GridIndices(std::shared_ptr< UserCtx > & user)
    {
      auto & info = user->info;
      xs = info.xs; xe = info.xs + info.xm;
      ys = info.ys; ye = info.ys + info.ym;
      zs = info.zs; ze = info.zs + info.zm;
      mx = info.mx; my = info.my; mz = info.mz;

      int i_begin = 1, i_end = mx-1;      // cross section in tecplot
      int j_begin = 1, j_end = my-1;
      int k_begin = 1, k_end = mz-1;
    
      PetscOptionsGetInt(PETSC_NULL, "-i_begin", &i_begin, PETSC_NULL);
      PetscOptionsGetInt(PETSC_NULL, "-i_end", &i_end, PETSC_NULL);
      PetscOptionsGetInt(PETSC_NULL, "-j_begin", &j_begin, PETSC_NULL);
      PetscOptionsGetInt(PETSC_NULL, "-j_end", &j_end, PETSC_NULL);
      PetscOptionsGetInt(PETSC_NULL, "-k_begin", &k_begin, PETSC_NULL);
      PetscOptionsGetInt(PETSC_NULL, "-k_end", &k_end, PETSC_NULL);

      xs = i_begin - 1, xe = i_end+1;
      ys = j_begin - 1, ye = j_end+1;
      zs = k_begin - 1, ze = k_end+1;

      IMax = i_end - i_begin + 1;
      JMax = j_end - j_begin + 1;
      KMax = k_end - k_begin + 1;
    }
};

PetscErrorCode TECInitialize(std::shared_ptr< UserCtx > & user,
                             char const * variables, char const * filename)
{
  INTEGER4 Debug = 0;
  INTEGER4 VIsDouble = 0;
  INTEGER4 DIsDouble = 0;
  INTEGER4 ZoneType = 0;
  INTEGER4 ICellMax = 0;
  INTEGER4 JCellMax = 0;
  INTEGER4 KCellMax = 0;
  INTEGER4 IsBlock = 1;
  INTEGER4 NumFaceConnections = 0;
  INTEGER4 FaceNeighborMode = 0;
  INTEGER4 ShareConnectivityFromZone = 0;
  // 1 is cell-centered 0 is node centered
  INTEGER4 LOC[40] = {1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  
  INTEGER4 I;
  auto idx = GridIndices(user);

  I = TECINI100((char*)"Flow", (char*)variables, (char*)filename, (char*)".",
                &Debug, &VIsDouble);
  if (I != 0)
  {
    printf("TECINI100 returned error %d\n", I);
    exit(1);
  }

  I = TECZNE100((char*)"Block 1",
                &ZoneType,      /* Ordered zone */
                &idx.IMax,
                &idx.JMax,
                &idx.KMax,
                &ICellMax,
                &JCellMax,
                &KCellMax,
                &IsBlock,       /* ISBLOCK  BLOCK format */
                &NumFaceConnections,
                &FaceNeighborMode,
                LOC,
                NULL,
                &ShareConnectivityFromZone); /* No connectivity sharing */
  if (I != 0)
  {
    printf("TECZNE100 returned error %d\n", I);
    exit(1);
  }
}

void TECOutputGrid(std::shared_ptr< UserCtx > & user)
{
  int i, j, k;
  GridIndices idx(user);
  INTEGER4 I, DIsDouble = 0, III = idx.IMax * idx.JMax * idx.KMax;;

  float * x = new float[III];

  for (k=idx.zs; k<idx.ze-1; k++)
  for (j=idx.ys; j<idx.ye-1; j++)
  for (i=idx.xs; i<idx.xe-1; i++)
    x[(k-idx.zs) * idx.IMax*idx.JMax + (j-idx.ys)*idx.IMax + (i-idx.xs)]
      = user->coor[k][j][i].x;
  I = TECDAT100(&III, &x[0], &DIsDouble);
  if (I != 0)
  {
    printf("TECDAT100 returned error %d\n", I);
    exit(1);
  }

  for (k=idx.zs; k<idx.ze-1; k++)
  for (j=idx.ys; j<idx.ye-1; j++)
  for (i=idx.xs; i<idx.xe-1; i++)
    x[(k-idx.zs) * idx.IMax*idx.JMax + (j-idx.ys)*idx.IMax + (i-idx.xs)]
      = user->coor[k][j][i].y;
  I = TECDAT100(&III, &x[0], &DIsDouble);
  if (I != 0)
  {
    printf("TECDAT100 returned error %d\n", I);
    exit(1);
  }

  for (k=idx.zs; k<idx.ze-1; k++)
  for (j=idx.ys; j<idx.ye-1; j++)
  for (i=idx.xs; i<idx.xe-1; i++)
    x[(k-idx.zs) * idx.IMax*idx.JMax + (j-idx.ys)*idx.IMax + (i-idx.xs)]
      = user->coor[k][j][i].z;
  I = TECDAT100(&III, &x[0], &DIsDouble);
  if (I != 0)
  {
    printf("TECDAT100 returned error %d\n", I);
    exit(1);
  }

  delete [] x;
}


void TECOutputVector(std::shared_ptr< UserCtx > & user, PetscReal ***vec)
{
  int i, j, k;
  GridIndices idx(user);
  INTEGER4 I, DIsDouble = 0, III = (idx.IMax-1)*(idx.JMax-1)*(idx.KMax-1);

  float * x = new float[III];

  for (k=idx.zs; k<idx.ze-2; k++)
  for (j=idx.ys; j<idx.ye-2; j++)
  for (i=idx.xs; i<idx.xe-2; i++)
    x[(k-idx.zs)*(idx.IMax-1)*(idx.JMax-1)+(j-idx.ys)*(idx.IMax-1)+(i-idx.xs)]
      = vec[k+1][j+1][i+1];

  I = TECDAT100(&III, &x[0], &DIsDouble);
  if (I != 0)
  {
    printf("TECDAT100 returned error %d\n", I);
    exit(1);
  }

  delete [] x;
}


// <<WRO>>
void TECOutputCmpnts(INTEGER4 IMax, INTEGER4 JMax, INTEGER4 KMax,
                     PetscInt xs, PetscInt xe, PetscInt ys, PetscInt ye,
                     PetscInt zs, PetscInt ze, Cmpnts ***xcat)
{
  PetscInt i, j, k;
  INTEGER4 I, DIsDouble=0;
  INTEGER4 III = (IMax-1)*(JMax-1)*(KMax-1);

  auto x = new float[III];

  if(!onlyV)
  {
    for (k=zs; k<ze-2; k++)
    for (j=ys; j<ye-2; j++)
    for (i=xs; i<xe-2; i++)
      x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)]
                                                       = xcat[k+1][j+1][i+1].x;
    I = TECDAT100(&III, &x[0], &DIsDouble);

    for (k=zs; k<ze-2; k++)
    for (j=ys; j<ye-2; j++)
    for (i=xs; i<xe-2; i++)
      x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)]
                                                       = xcat[k+1][j+1][i+1].y;
    I = TECDAT100(&III, &x[0], &DIsDouble);

    for (k=zs; k<ze-2; k++)
    for (j=ys; j<ye-2; j++)
    for (i=xs; i<xe-2; i++)
      x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)]
                                                       = xcat[k+1][j+1][i+1].z;
    I = TECDAT100(&III, &x[0], &DIsDouble);
  }

  for (k=zs; k<ze-2; k++)
  for (j=ys; j<ye-2; j++)
  for (i=xs; i<xe-2; i++)
    x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] =
      sqrt(xcat[k+1][j+1][i+1].x*xcat[k+1][j+1][i+1].x
         + xcat[k+1][j+1][i+1].y*xcat[k+1][j+1][i+1].y
         + xcat[k+1][j+1][i+1].z*xcat[k+1][j+1][i+1].z);
  I = TECDAT100(&III, &x[0], &DIsDouble);

  delete [] x;
}

void calcCellCenterValuesCmpnts(int mx, int my, int mz,
                                Cmpnts ***dst, Cmpnts  ***src)
{
  PetscInt i, j, k;

  for (k=0; k<mz-1; k++)
  for (j=0; j<my-1; j++)
  for (i=0; i<mx-1; i++)
  {
    dst[k][j][i].x = 0.125 *
      (  src[k  ][j  ][i].x + src[k  ][j  ][i+1].x
       + src[k  ][j+1][i].x + src[k  ][j+1][i+1].x
       + src[k+1][j  ][i].x + src[k+1][j  ][i+1].x
       + src[k+1][j+1][i].x + src[k+1][j+1][i+1].x);
    dst[k][j][i].y = 0.125 *
      (  src[k  ][j  ][i].y + src[k  ][j  ][i+1].y
       + src[k  ][j+1][i].y + src[k  ][j+1][i+1].y
       + src[k+1][j  ][i].y + src[k+1][j  ][i+1].y
       + src[k+1][j+1][i].y + src[k+1][j+1][i+1].y);
    dst[k][j][i].z = 0.125 *
      (  src[k  ][j  ][i].z + src[k  ][j  ][i+1].z
       + src[k  ][j+1][i].z + src[k  ][j+1][i+1].z
       + src[k+1][j  ][i].z + src[k+1][j  ][i+1].z
       + src[k+1][j+1][i].z + src[k+1][j+1][i+1].z);
  }
}


PetscErrorCode TECIOOut_V(std::shared_ptr< UserCtx > & user)
{
  PetscInt bi;

  char filen[80];
  sprintf(filen, "%sResult%06d.plt", prefix, ti);

  printf("\nGenerating file %s\n", filen);

  // WRO - paraview
  if (paraview)
  {
    char paran[strlen(prefix) + 20];
    sprintf(paran, "%sparaview%06d.vtk", prefix, ti);
    Paraview::init(paran);
  }

  INTEGER4 I, Debug, VIsDouble, DIsDouble, III, IMax, JMax, KMax;
  VIsDouble = 0;
  DIsDouble = 0;
  Debug = 0;



  if(onlyV)   {
    if(cs) I = TECINI100((char*)"Flow", (char*)"X Y Z UU Cs", filen, (char*)".", &Debug, &VIsDouble);
    else if(onlyV==2) I = TECINI100((char*)"Flow", (char*)"X Y Z UU", filen, (char*)".", &Debug, &VIsDouble);
    else I = TECINI100((char*)"Flow", (char*)"X Y Z UU Nv", filen, (char*)".", &Debug, &VIsDouble);
  }
  else if(rans/* && rans_output*/) {
    if(conv_diff) I = TECINI100((char*)"Flow", (char*)"X Y Z U V W UU P Nv K Omega Nut C ", filen, (char*)".", &Debug, &VIsDouble);
    else if(levelset) I = TECINI100((char*)"Flow", (char*)"X Y Z U V W UU P Nv K Omega Nut Level Fr", filen, (char*)".", &Debug, &VIsDouble);
    else I = TECINI100((char*)"Flow", (char*)"X Y Z U V W UU P Nv K Omega Nut ", filen, (char*)".", &Debug, &VIsDouble);
  }
  else  {
    if(cs) I = TECINI100((char*)"Flow", (char*)"X Y Z U V W UU P Cs", filen, (char*)".", &Debug, &VIsDouble);
    else if (levelset && !conv_diff) I = TECINI100((char*)"Flow", (char*)"X Y Z U V W UU P Nv Level Fr", filen, (char*)".", &Debug, &VIsDouble);
    else if (conv_diff && !sandwave && !levelset)I = TECINI100((char*)"Flow", (char*)"X Y Z U V W UU P Nv C", filen, (char*)".", &Debug, &VIsDouble);
    else if (sandwave && conv_diff)I = TECINI100((char*)"Flow", (char*)"X Y Z U V W UU P Nv ustar C", filen, (char*)".", &Debug, &VIsDouble);
    else if (conv_diff && levelset)I = TECINI100((char*)"Flow", (char*)"X Y Z U V W UU P Nv C Level Fr", filen, (char*)".", &Debug, &VIsDouble);
    else I = TECINI100((char*)"Flow", (char*)"X Y Z U V W UU P Nv", filen, (char*)".", &Debug, &VIsDouble);

  }

  for (bi=0; bi<block_number; bi++) {
    DA da = /* user[bi]. */ user->da, fda = /* user[bi]. */ user->fda;
    DALocalInfo info = /* user[bi]. */ user->info;

    PetscInt        xs = info.xs, xe = info.xs + info.xm;
    PetscInt        ys = info.ys, ye = info.ys + info.ym;
    PetscInt        zs = info.zs, ze = info.zs + info.zm;
    PetscInt        mx = info.mx, my = info.my, mz = info.mz;

    PetscInt        lxs, lys, lzs, lxe, lye, lze;
    PetscInt        i, j, k;
    PetscReal       ***aj;
    Cmpnts  ***ucat, ***coor, ***ucat_o, ***csi, ***eta, ***zet;
    PetscReal       ***p, ***nvert, ***level;
    Vec             Coor, zCoor, nCoor;
    VecScatter      ctx;
    Vec K_Omega;

    DAVecGetArray(/* user[bi]. */ user->da, /* user[bi]. */ user->Nvert, &nvert);
    // DAVecGetArray(/* user[bi]. */ user->da, /* user[bi]. */ user->Aj, &aj);
    // DAVecGetArray(/* user[bi]. */ user->fda, /* user[bi]. */ user->Csi, &csi);
    // DAVecGetArray(/* user[bi]. */ user->fda, /* user[bi]. */ user->Eta, &eta);
    // DAVecGetArray(/* user[bi]. */ user->fda, /* user[bi]. */ user->Zet, &zet);
    DAVecGetArray(/* user[bi]. */ user->fda, /* user[bi]. */ user->Ucat, &ucat);

    INTEGER4        ZoneType=0, ICellMax=0, JCellMax=0, KCellMax=0;
    INTEGER4        IsBlock=1, NumFaceConnections=0, FaceNeighborMode=0;
    INTEGER4    ShareConnectivityFromZone=0;
    INTEGER4        LOC[40] = {1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}; /* 1 is cell-centered 0 is node centered */


    /*************************/
    printf("mi=%d, mj=%d, mk=%d\n", mx, my, mz);
    printf("xs=%d, xe=%d\n", xs, xe);
    printf("ys=%d, ye=%d\n", ys, ye);
    printf("zs=%d, ze=%d\n", zs, ze);
    //exit(0);

    i_begin = 1, i_end = mx-1;      // cross section in tecplot
    j_begin = 1, j_end = my-1;
    k_begin = 1, k_end = mz-1;

    PetscOptionsGetInt(PETSC_NULL, "-i_begin", &i_begin, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-i_end", &i_end, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-j_begin", &j_begin, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-j_end", &j_end, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-k_begin", &k_begin, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, "-k_end", &k_end, PETSC_NULL);

    xs = i_begin - 1, xe = i_end+1;
    ys = j_begin - 1, ye = j_end+1;
    zs = k_begin - 1, ze = k_end+1;

    printf("xs=%d, xe=%d\n", xs, xe);
    printf("ys=%d, ye=%d\n", ys, ye);
    printf("zs=%d, ze=%d\n", zs, ze);
    //exit(0);
    //xs=0, xe=nsection+1;
    /*************************/


    // WRO - output max values of vmag
    {
      StoreMax<> maxs("vmag");
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++)
          {
            double x = sqrt(ucat[k+1][j+1][i+1].x*ucat[k+1][j+1][i+1].x+ucat[k+1][j+1][i+1].y*ucat[k+1][j+1][i+1].y+ucat[k+1][j+1][i+1].z*ucat[k+1][j+1][i+1].z);
            maxs.checkEntry(i, j, k, x);
          }
      maxs.output();
    }


    if (vc) {
      LOC[3]=0; LOC[4]=0; LOC[5]=0; LOC[6]=0;
    }
    else if(onlyV) {
      LOC[4]=0; LOC[5]=0; LOC[6]=0;
    }
    /*
      IMax = mx-1;
      JMax = my-1;
      KMax = mz-1;
    */

    IMax = i_end - i_begin + 1;
    JMax = j_end - j_begin + 1;
    KMax = k_end - k_begin + 1;

    I = TECZNE100((char*)"Block 1",
                  &ZoneType,      /* Ordered zone */
                  &IMax,
                  &JMax,
                  &KMax,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &IsBlock,       /* ISBLOCK  BLOCK format */
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  LOC,
                  NULL,
                  &ShareConnectivityFromZone); /* No connectivity sharing */

    //III = (mx-1) * (my-1) * (mz-1);
    III = IMax*JMax*KMax;

    DAGetCoordinates(da, &Coor);
    DAVecGetArray(fda, Coor, &coor);

    float *x;
    x = new float [III];

    for (k=zs; k<ze-1; k++)
      for (j=ys; j<ye-1; j++)
        for (i=xs; i<xe-1; i++) {
          x[(k-zs) * IMax*JMax + (j-ys)*IMax + (i-xs)] = coor[k][j][i].x;
        }
    I = TECDAT100(&III, &x[0], &DIsDouble);

    for (k=zs; k<ze-1; k++)
      for (j=ys; j<ye-1; j++)
        for (i=xs; i<xe-1; i++) {
          x[(k-zs) * IMax*JMax + (j-ys)*IMax + (i-xs)] = coor[k][j][i].y;
        }
    I = TECDAT100(&III, &x[0], &DIsDouble);

    for (k=zs; k<ze-1; k++)
      for (j=ys; j<ye-1; j++)
        for (i=xs; i<xe-1; i++) {
          x[(k-zs) * IMax*JMax + (j-ys)*IMax + (i-xs)] = coor[k][j][i].z;
        }
    I = TECDAT100(&III, &x[0], &DIsDouble);


    // WRO - output data in vtk format for paraview.
    if (paraview)
    {
      printf("writing coor in paraview format\n");
      Paraview::ptr()->outputGrid(xs, xe, ys, ye, zs, ze, coor);
    }


    DAVecRestoreArray(fda, Coor, &coor);
    delete []x;


    if(!vc) {
      DAVecGetArray(/* user[bi]. */ user->fda,
                    /* user[bi]. */ user->Ucat_o, &ucat_o);
      
      calcCellCenterValuesCmpnts(mx, my, mz, ucat_o, ucat);
      
      for (k=zs; k<ze-1; k++)
        for (j=ys; j<ye-1; j++)
          for (i=xs; i<xe-1; i++) {
            x[k * (mx-1)*(my-1) + j*(mx-1) + i] = ucat_o[k][j][i].x;
          }
      if(!onlyV) I = TECDAT100(&III, &x[0], &DIsDouble);

      for (k=zs; k<ze-1; k++)
        for (j=ys; j<ye-1; j++)
          for (i=xs; i<xe-1; i++) {
            x[k * (mx-1)*(my-1) + j*(mx-1) + i] = ucat_o[k][j][i].y;
          }
      if(!onlyV) I = TECDAT100(&III, &x[0], &DIsDouble);

      for (k=zs; k<ze-1; k++)
        for (j=ys; j<ye-1; j++)
          for (i=xs; i<xe-1; i++) {
            x[k * (mx-1)*(my-1) + j*(mx-1) + i] = ucat_o[k][j][i].z;
          }

      if(!onlyV) I = TECDAT100(&III, &x[0], &DIsDouble);

      for (k=zs; k<ze-1; k++)
        for (j=ys; j<ye-1; j++)
          for (i=xs; i<xe-1; i++) {
            x[k * (mx-1)*(my-1) + j*(mx-1) + i] = sqrt( ucat_o[k][j][i].x*ucat_o[k][j][i].x + ucat_o[k][j][i].y*ucat_o[k][j][i].y + ucat_o[k][j][i].z*ucat_o[k][j][i].z );
          }
      I = TECDAT100(&III, &x[0], &DIsDouble);

      DAVecRestoreArray(/* user[bi]. */ user->fda, /* user[bi]. */ user->Ucat_o, &ucat_o);
      delete []x;
    }
    else
      // <<WRO>>
      TECOutputCmpnts(IMax, JMax, KMax, xs, xe, ys, ye, zs, ze, ucat);


    III = (IMax-1)*(JMax-1)*(KMax-1);
    //III = (mx-2) * (my-2) * (mz-2);
    x = new float [III];
    //x.resize (III);

    if(!onlyV) {
      DAVecGetArray(/* user[bi]. */ user->da, /* user[bi]. */ user->P, &p);
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) {
            x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = p[k+1][j+1][i+1];
          }
      I = TECDAT100(&III, &x[0], &DIsDouble);
      DAVecRestoreArray(/* user[bi]. */ user->da, /* user[bi]. */ user->P, &p);
    }


    if(onlyV!=2) {
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) {
            x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = nvert[k+1][j+1][i+1];
          }
      I = TECDAT100(&III, &x[0], &DIsDouble);
    }

    if(onlyV==2) { // Z Vorticity
      /*
        for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
        for (i=xs; i<xe-2; i++) {
        double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
        double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
        double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
        double eta0= eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
        double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
        double ajc = aj[k][j][i];

        if(i==0 || j==0 || k==0) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = 0;
        else {
        Compute_du_center (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
        Compute_du_dxyz (       csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
        &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
        x[k * (mx-2)*(my-2) + j*(mx-2) + i] = dv_dx - du_dy;
        }
        }
        I = TECDAT100(&III, &x[0], &DIsDouble);
      */
    }

    if(!onlyV && rans /*&& rans_output*/) {
      Cmpnts2 ***komega;
      DACreateGlobalVector(user->fda2, &K_Omega);
      PetscViewer     viewer;
      sprintf(filen, "kfield%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoadIntoVector(viewer, K_Omega);
      PetscViewerDestroy(viewer);
      DAVecGetArray(/* user[bi]. */ user->fda2, K_Omega, &komega);

      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = komega[k+1][j+1][i+1].x;
      I = TECDAT100(&III, &x[0], &DIsDouble);

      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = komega[k+1][j+1][i+1].y;
      I = TECDAT100(&III, &x[0], &DIsDouble);

      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = komega[k+1][j+1][i+1].x/(komega[k+1][j+1][i+1].y+1.e-20);
      I = TECDAT100(&III, &x[0], &DIsDouble);

      DAVecRestoreArray(/* user[bi]. */ user->fda2, K_Omega, &komega);
      VecDestroy(K_Omega);
    }

    if(conv_diff && sandwave) {
      PetscReal ***ustar_;
      Vec lUstar_;
      DACreateGlobalVector(user->da, &lUstar_);
      PetscViewer     viewer;
      sprintf(filen, "ustarfield%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoadIntoVector(viewer, lUstar_);
      PetscViewerDestroy(viewer);
      DAVecGetArray(/* user[bi]. */ user->da, lUstar_, &ustar_);

      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = ustar_[k+1][j+1][i+1];

      if(width_ave_instant) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz); // Added to compute the Hacker Case for the width averaged density current June 2012

      I = TECDAT100(&III, &x[0], &DIsDouble);

      DAVecRestoreArray(/* user[bi]. */ user->da, lUstar_, &ustar_);
      VecDestroy(lUstar_);
    }


    if(conv_diff) {
      PetscReal ***conc;
      Vec Conc;
      DACreateGlobalVector(user->da, &Conc);
      PetscViewer     viewer;
      sprintf(filen, "cfield%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoadIntoVector(viewer, Conc);
      PetscViewerDestroy(viewer);
      DAVecGetArray(/* user[bi]. */ user->da, Conc, &conc);

      {
        StoreMax<20> maxc("conc");
        for (k=zs; k<ze-2; k++)
          for (j=ys; j<ye-2; j++)
            for (i=xs; i<xe-2; i++)
            {
              maxc.checkEntry(i, j, k, conc[k+1][j+1][i+1]);
              x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = conc[k+1][j+1][i+1];
            }
        maxc.output();
      }

      if(width_ave_instant) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz); // Added to compute the Hacker Case for the width averaged density current June 2012

      I = TECDAT100(&III, &x[0], &DIsDouble);


      // WRO - write concentration data to paraview file.
      if (paraview)
      {
        printf("writing concentration to paraview format\n");
        Paraview::ptr()->outputData("C", &x[0], III);
      }


      DAVecRestoreArray(/* user[bi]. */ user->da, Conc, &conc);
      VecDestroy(Conc);
    }

    if(levelset) {
      PetscReal ***level;
      Vec Levelset;
      DACreateGlobalVector(user->da, &Levelset);
      PetscViewer     viewer;
      sprintf(filen, "lfield%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoadIntoVector(viewer, Levelset);
      PetscViewerDestroy(viewer);
      DAVecGetArray(/* user[bi]. */ user->da, Levelset, &level);

      // level
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = level[k+1][j+1][i+1];
      I = TECDAT100(&III, &x[0], &DIsDouble);

      // Froude number
      for (k=zs; k<ze-2; k++)
        for (i=xs; i<xe-2; i++) {
          double Fr=0; // Froude number
          double G=0;
          double min_phi=1.e10, max_phi=-1.e10;
          double sum_velocity_magnitude = 0;
          int count = 0;

          if( fabs(gravity_x)>1.e-10 ) G=fabs(gravity_x);
          else if( fabs(gravity_y)>1.e-10 ) G=fabs(gravity_y);
          else if( fabs(gravity_z)>1.e-10 ) G=fabs(gravity_z);

          for (j=ys; j<ye-2; j++) {
            if(level[k+1][j+1][i+1]>0 && nvert[k+1][j+1][i+1]<0.1) {
              min_phi = std::min (min_phi, level[k+1][j+1][i+1]);
              max_phi = std::max (max_phi, level[k+1][j+1][i+1]);
              sum_velocity_magnitude += sqrt(ucat[k+1][j+1][i+1].x*ucat[k+1][j+1][i+1].x+ucat[k+1][j+1][i+1].y*ucat[k+1][j+1][i+1].y+ucat[k+1][j+1][i+1].z*ucat[k+1][j+1][i+1].z);
              count += 1;
            }
          }

          if(count) {
            double UU = sum_velocity_magnitude / (double)count;
            double depth = max_phi - min_phi;
            Fr = UU / sqrt ( G * depth );
            if(count<=2 || Fr>1000) Fr=0;
          }

          for (j=ys; j<ye-2; j++) x[ (k-zs) * (IMax-1)*(JMax-1) + (j-ys) * (IMax-1) + (i-xs)] = Fr;
        }
      I = TECDAT100(&III, &x[0], &DIsDouble);
      DAVecRestoreArray(/* user[bi]. */ user->da, Levelset, &level);
      VecDestroy(Levelset);
    }

    delete []x;

    // DAVecRestoreArray(/* user[bi]. */ user->da, /* user[bi]. */ user->Aj, &aj);
    // DAVecRestoreArray(/* user[bi]. */ user->fda, /* user[bi]. */ user->Csi, &csi);
    // DAVecRestoreArray(/* user[bi]. */ user->fda, /* user[bi]. */ user->Eta, &eta);
    // DAVecRestoreArray(/* user[bi]. */ user->fda, /* user[bi]. */ user->Zet, &zet);
    DAVecRestoreArray(/* user[bi]. */ user->fda, /* user[bi]. */ user->Ucat, &ucat);
    DAVecRestoreArray(/* user[bi]. */ user->da, /* user[bi]. */ user->Nvert, &nvert);
  }
  I = TECEND100();
  return 0;
}


PetscErrorCode TECIOOut_Averaging(std::shared_ptr< UserCtx > & user)
{
  PetscInt bi;

  char filen[80];
  sprintf(filen, "%sResult%06d-avg.plt", prefix, ti);

  INTEGER4 I, Debug, VIsDouble, DIsDouble, III, IMax, JMax, KMax;
  VIsDouble = 0;
  DIsDouble = 0;
  Debug = 0;

  if(pcr) I = TECINI100((char*)"Averaging", (char*)"X Y Z P Velocity_Magnitude Nv",  filen, (char*)".",  &Debug,  &VIsDouble);
  else if(avg==1) {

    if(averaging_option==1) {
      //I = TECINI100((char*)"Averaging", "X Y Z U V W uu vv ww uv vw uw UV_ VW_ UW_ Nv",  filen, (char*)".",  &Debug,  &VIsDouble); //OSL
      //I = TECINI100((char*)"Averaging", "X Y Z U V W uu vv ww uv vw uw UU_ VV_ WW_ Nv",  filen, (char*)".",  &Debug,  &VIsDouble); //OSL
      if(levelset)I = TECINI100((char*)"Averaging", (char*)"X Y Z U V W uu vv ww uv vw uw l lrms K Nv",  filen, (char*)".",  &Debug,  &VIsDouble);
      if(!levelset)I = TECINI100((char*)"Averaging", (char*)"X Y Z U V W uu vv ww uv vw uw K Nv",  filen, (char*)".",  &Debug,  &VIsDouble);

    }
    else if(averaging_option==2) {
      //I = TECINI100((char*)"Averaging", "X Y Z U V W uu vv ww uv vw uw  UV_ VW_ UW_ P pp Nv",  filen, (char*)".",  &Debug,  &VIsDouble);
      I = TECINI100((char*)"Averaging", (char*)"X Y Z U V W uu vv ww uv vw uw  P pp Nv",  filen, (char*)".",  &Debug,  &VIsDouble);
    }
    else if(averaging_option==3) {
      //I = TECINI100((char*)"Averaging", "X Y Z U V W uu vv ww uv vw uw UV_ VW_ UW_ P pp Vortx Vorty Vortz vortx2 vorty2 vortz2 Nv",  filen, (char*)".",  &Debug,  &VIsDouble);
      I = TECINI100((char*)"Averaging", (char*)"X Y Z U V W uu vv ww uv vw uw P pp Vortx Vorty Vortz vortx2 vorty2 vortz2 Nv",  filen, (char*)".",  &Debug,  &VIsDouble);
    }

  }
  else if(avg==2 && !sandwave){
    if(conv_diff && levelset) I = TECINI100((char*)"Averaging", (char*)"X Y Z U V W K C l lrms Nv",  filen, (char*)".",  &Debug,  &VIsDouble);
    else if(conv_diff) I = TECINI100((char*)"Averaging", (char*)"X Y Z U V W K C Nv",  filen, (char*)".",  &Debug,  &VIsDouble);
    else if(levelset) I = TECINI100((char*)"Averaging", (char*)"X Y Z U V W K l lrms Nv",  filen, (char*)".",  &Debug,  &VIsDouble);
    else I = TECINI100((char*)"Averaging", (char*)"X Y Z U V W K Nv",  filen, (char*)".",  &Debug,  &VIsDouble);
  }
  else if(avg==2 && sandwave){
    I = TECINI100((char*)"Averaging", (char*)"X Y Z uf vf wf F ustar ustarf Nv",  filen, (char*)".",  &Debug,  &VIsDouble);
  }

  for (bi=0; bi<block_number; bi++) {
    DA da = /* user[bi]. */ user->da, fda = /* user[bi]. */ user->fda;
    DALocalInfo info = /* user[bi]. */ user->info;

    PetscInt        xs = info.xs, xe = info.xs + info.xm;
    PetscInt        ys = info.ys, ye = info.ys + info.ym;
    PetscInt        zs = info.zs, ze = info.zs + info.zm;
    PetscInt        mx = info.mx, my = info.my, mz = info.mz;

    PetscInt        lxs, lys, lzs, lxe, lye, lze;
    PetscInt        i, j, k;
    PetscReal       ***aj;
    Cmpnts  ***ucat, ***coor, ***ucat_o;
    Cmpnts  ***u2sum, ***u1sum,  ***usum;
    PetscReal       ***p, ***nvert;
    Vec             Coor, zCoor, nCoor;

    //VecScatter ctx;

    Vec X, Y, Z, U, V, W;

    INTEGER4        ZoneType=0, ICellMax=0, JCellMax=0, KCellMax=0;
    INTEGER4        IsBlock=1, NumFaceConnections=0, FaceNeighborMode=0;
    INTEGER4    ShareConnectivityFromZone=0;
    INTEGER4        LOC[100] = {1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    /* 1 is cell-centered   0 is node centered */

    IMax = mx-1; JMax = my-1; KMax = mz-1;

    I = TECZNE100((char*)"Block 1",
                  &ZoneType,      /* Ordered zone */
                  &IMax,
                  &JMax,
                  &KMax,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &IsBlock,       /* ISBLOCK  BLOCK format */
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  LOC,
                  NULL,
                  &ShareConnectivityFromZone); /* No connectivity sharing */

    float *x;

    x = new float [(mx-1)*(my-1)*(mz-1)];
    III = (mx-1) * (my-1) * (mz-1);

    DAGetCoordinates(da, &Coor);
    DAVecGetArray(fda, Coor, &coor);

    // X
    for (k=zs; k<ze-1; k++)
      for (j=ys; j<ye-1; j++)
        for (i=xs; i<xe-1; i++) x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].x;
    I = TECDAT100(&III, x, &DIsDouble);

    // Y
    for (k=zs; k<ze-1; k++)
      for (j=ys; j<ye-1; j++)
        for (i=xs; i<xe-1; i++) x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].y;
    I = TECDAT100(&III, x, &DIsDouble);

    // Z
    for (k=zs; k<ze-1; k++)
      for (j=ys; j<ye-1; j++)
        for (i=xs; i<xe-1; i++) x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].z;
    I = TECDAT100(&III, x, &DIsDouble);

    DAVecRestoreArray(fda, Coor, &coor);

    //delete []x;
    double N=(double)ti+1.0;
    //x = new float [(mx-2)*(my-2)*(mz-2)];

    III = (mx-2) * (my-2) * (mz-2);

    if(pcr)  {
      DAVecGetArray(/* user[bi]. */ user->da, /* user[bi]. */ user->P, &p);
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = p[k+1][j+1][i+1];
      I = TECDAT100(&III, x, &DIsDouble);
      DAVecRestoreArray(/* user[bi]. */ user->da, /* user[bi]. */ user->P, &p);

      // Load ucat
      PetscViewer     viewer;
      sprintf(filen, "ufield%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoadIntoVector(viewer, (/* user[bi]. */ user->Ucat));
      PetscViewerDestroy(viewer);

      DAVecGetArray(/* user[bi]. */ user->fda, /* user[bi]. */ user->Ucat, &ucat);
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i]
                                    = sqrt ( ucat[k+1][j+1][i+1].x*ucat[k+1][j+1][i+1].x + ucat[k+1][j+1][i+1].y*ucat[k+1][j+1][i+1].y + ucat[k+1][j+1][i+1].z*ucat[k+1][j+1][i+1].z );
      I = TECDAT100(&III, x, &DIsDouble);
      DAVecRestoreArray(/* user[bi]. */ user->fda, /* user[bi]. */ user->Ucat, &ucat);
    }
    else if(avg==1) {

      Vec K_sum;
      PetscReal ***ksum;
      Vec Level_sum;
      Vec Level_square_sum;
      PetscReal ***level_sum;
      PetscReal ***level2sum;

      PetscViewer viewer;
      char filen[128];


      DACreateGlobalVector(/* user[bi]. */ user->fda, &/* user[bi]. */ user->Ucat_sum);
      DACreateGlobalVector(/* user[bi]. */ user->fda, &/* user[bi]. */ user->Ucat_cross_sum);
      DACreateGlobalVector(/* user[bi]. */ user->fda, &/* user[bi]. */ user->Ucat_square_sum);
      if(levelset) {
        DACreateGlobalVector(user->da, &Level_sum);
        DACreateGlobalVector(user->da, &Level_square_sum);}

      if(levelset) {
        sprintf(filen, "slevel_%06d_%1d.dat", ti, user->_this);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
        VecLoadIntoVector(viewer, Level_sum);
        PetscViewerDestroy(viewer);

        sprintf(filen, "slevel2_%06d_%1d.dat", ti, user->_this);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
        VecLoadIntoVector(viewer, Level_square_sum);
        PetscViewerDestroy(viewer);
      }

      sprintf(filen, "su0_%06d_%1d.dat", ti, /* user[bi]. */ user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoadIntoVector(viewer, (/* user[bi]. */ user->Ucat_sum));
      PetscViewerDestroy(viewer);

      sprintf(filen, "su1_%06d_%1d.dat", ti, /* user[bi]. */ user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoadIntoVector(viewer, (/* user[bi]. */ user->Ucat_cross_sum));
      PetscViewerDestroy(viewer);

      sprintf(filen, "su2_%06d_%1d.dat", ti, /* user[bi]. */ user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoadIntoVector(viewer, (/* user[bi]. */ user->Ucat_square_sum));
      PetscViewerDestroy(viewer);

      DAVecGetArray(/* user[bi]. */ user->fda, /* user[bi]. */ user->Ucat_sum, &usum);
      DAVecGetArray(/* user[bi]. */ user->fda, /* user[bi]. */ user->Ucat_cross_sum, &u1sum);
      DAVecGetArray(/* user[bi]. */ user->fda, /* user[bi]. */ user->Ucat_square_sum, &u2sum);
      if(levelset) {
        DAVecGetArray(/* user[bi]. */ user->da, Level_sum, &level_sum);
        DAVecGetArray(/* user[bi]. */ user->da, Level_square_sum, &level2sum);
      }

      // U
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) {x[k * (mx-2)*(my-2) + j*(mx-2) + i] = usum[k+1][j+1][i+1].x/N;}
      if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      I = TECDAT100(&III, x, &DIsDouble);

      // V
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++){x[k * (mx-2)*(my-2) + j*(mx-2) + i] = usum[k+1][j+1][i+1].y/N;}
      if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      I = TECDAT100(&III, x, &DIsDouble);

      // W
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) {x[k * (mx-2)*(my-2) + j*(mx-2) + i] = usum[k+1][j+1][i+1].z/N;}
      if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      I = TECDAT100(&III, x, &DIsDouble);

      // uu, u rms
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for(i=xs; i<xe-2; i++) {double U = usum[k+1][j+1][i+1].x/N;
            x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( u2sum[k+1][j+1][i+1].x/N - U*U );
          }
      if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      I = TECDAT100(&III, x, &DIsDouble);

      // vv, v rms
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) {
            double V = usum[k+1][j+1][i+1].y/N;
            x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( u2sum[k+1][j+1][i+1].y/N - V*V );
          }
      if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      I = TECDAT100(&III, x, &DIsDouble);

      // ww, w rms
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) {
            double W = usum[k+1][j+1][i+1].z/N;
            x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( u2sum[k+1][j+1][i+1].z/N - W*W );
          }
      if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      I = TECDAT100(&III, x, &DIsDouble);

      // uv
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) {
            double UV = usum[k+1][j+1][i+1].x*usum[k+1][j+1][i+1].y / (N*N);
            x[k * (mx-2)*(my-2) + j*(mx-2) + i] = u1sum[k+1][j+1][i+1].x/N - UV;
          }
      if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      I = TECDAT100(&III, x, &DIsDouble);

      // vw
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) {
            double VW = usum[k+1][j+1][i+1].y*usum[k+1][j+1][i+1].z / (N*N);
            x[k * (mx-2)*(my-2) + j*(mx-2) + i] = u1sum[k+1][j+1][i+1].y/N - VW;
          }
      if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      I = TECDAT100(&III, x, &DIsDouble);

      // wu
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) {
            double WU = usum[k+1][j+1][i+1].z*usum[k+1][j+1][i+1].x / (N*N);
            x[k * (mx-2)*(my-2) + j*(mx-2) + i] = u1sum[k+1][j+1][i+1].z/N - WU;
          }
      if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      I = TECDAT100(&III, x, &DIsDouble);

      // Level
      if(levelset){
        for (k=zs; k<ze-2; k++)
          for (j=ys; j<ye-2; j++)
            for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = level_sum[k+1][j+1][i+1]/N; // l
        if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        I = TECDAT100(&III, x, &DIsDouble);

        for (k=zs; k<ze-2; k++)
          for (j=ys; j<ye-2; j++)
            for (i=xs; i<xe-2; i++){
              double L = level_sum[k+1][j+1][i+1]/N;
              x[k * (mx-2)*(my-2) + j*(mx-2) + i] = sqrt( level2sum[k+1][j+1][i+1]/N - L*L ); // l rms
            }
        if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        I = TECDAT100(&III, x, &DIsDouble);
      }

      // k
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for(i=xs; i<xe-2; i++) {
            if(rans)  {
              x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ksum[k+1][j+1][i+1]/N;
            }
            else {
              double U = usum[k+1][j+1][i+1].x/N;
              double V = usum[k+1][j+1][i+1].y/N;
              double W = usum[k+1][j+1][i+1].z/N;
              x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( u2sum[k+1][j+1][i+1].x/N - U*U ) + ( u2sum[k+1][j+1][i+1].y/N - V*V ) + ( u2sum[k+1][j+1][i+1].z/N - W*W );
              x[k * (mx-2)*(my-2) + j*(mx-2) + i] *= 0.5;
            }
          }
      if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
      I = TECDAT100(&III, x, &DIsDouble);

      if(levelset){
        DAVecRestoreArray(/* user[bi]. */ user->da, Level_sum, &level_sum);
        DAVecRestoreArray(/* user[bi]. */ user->da, Level_square_sum, &level2sum);
      }

      if(levelset) {
        VecDestroy(Level_sum);
        VecDestroy(Level_square_sum);
      }

      if(averaging_option>=2) {
        Vec P_sum, P_square_sum;
        PetscReal ***psum, ***p2sum;

        DACreateGlobalVector(/* user[bi]. */ user->da, &P_sum);
        DACreateGlobalVector(/* user[bi]. */ user->da, &P_square_sum);

        sprintf(filen, "sp_%06d_%1d.dat", ti, /* user[bi]. */ user->_this);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
        VecLoadIntoVector(viewer, P_sum);
        PetscViewerDestroy(viewer);

        sprintf(filen, "sp2_%06d_%1d.dat", ti, /* user[bi]. */ user->_this);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
        VecLoadIntoVector(viewer, P_square_sum);
        PetscViewerDestroy(viewer);

        DAVecGetArray(/* user[bi]. */ user->da, P_sum, &psum);
        DAVecGetArray(/* user[bi]. */ user->da, P_square_sum, &p2sum);

        for (k=zs; k<ze-2; k++)
          for (j=ys; j<ye-2; j++)
            for (i=xs; i<xe-2; i++) {
              double P = psum[k+1][j+1][i+1]/N;
              x[k * (mx-2)*(my-2) + j*(mx-2) + i] = P;
            }
        if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        I = TECDAT100(&III, x, &DIsDouble);


        for (k=zs; k<ze-2; k++)
          for (j=ys; j<ye-2; j++)
            for (i=xs; i<xe-2; i++) {
              double P = psum[k+1][j+1][i+1]/N;
              x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( p2sum[k+1][j+1][i+1]/N - P*P );
            }
        if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        I = TECDAT100(&III, x, &DIsDouble);

        DAVecRestoreArray(/* user[bi]. */ user->da, P_sum, &psum);
        DAVecRestoreArray(/* user[bi]. */ user->da, P_square_sum, &p2sum);

        VecDestroy(P_sum);
        VecDestroy(P_square_sum);
      }

      if(averaging_option>=3) {
        Vec Vort_sum, Vort_square_sum;
        Cmpnts ***vortsum, ***vort2sum;

        DACreateGlobalVector(/* user[bi]. */ user->fda, &Vort_sum);
        DACreateGlobalVector(/* user[bi]. */ user->fda, &Vort_square_sum);

        sprintf(filen, "svo_%06d_%1d.dat", ti, /* user[bi]. */ user->_this);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
        VecLoadIntoVector(viewer, Vort_sum);
        PetscViewerDestroy(viewer);

        sprintf(filen, "svo2_%06d_%1d.dat", ti, /* user[bi]. */ user->_this);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
        VecLoadIntoVector(viewer, Vort_square_sum);
        PetscViewerDestroy(viewer);

        DAVecGetArray(/* user[bi]. */ user->fda, Vort_sum, &vortsum);
        DAVecGetArray(/* user[bi]. */ user->fda, Vort_square_sum, &vort2sum);

        for (k=zs; k<ze-2; k++)
          for (j=ys; j<ye-2; j++)
            for (i=xs; i<xe-2; i++) {
              double vortx = vortsum[k+1][j+1][i+1].x/N;
              x[k * (mx-2)*(my-2) + j*(mx-2) + i] = vortx;
            }
        if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        I = TECDAT100(&III, x, &DIsDouble);

        for (k=zs; k<ze-2; k++)
          for (j=ys; j<ye-2; j++)
            for (i=xs; i<xe-2; i++) {
              double vorty = vortsum[k+1][j+1][i+1].y/N;
              x[k * (mx-2)*(my-2) + j*(mx-2) + i] = vorty;
            }
        if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        I = TECDAT100(&III, x, &DIsDouble);

        for (k=zs; k<ze-2; k++)
          for (j=ys; j<ye-2; j++)
            for (i=xs; i<xe-2; i++) {
              double vortz = vortsum[k+1][j+1][i+1].z/N;
              x[k * (mx-2)*(my-2) + j*(mx-2) + i] = vortz;
            }
        if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        I = TECDAT100(&III, x, &DIsDouble);

        for (k=zs; k<ze-2; k++)
          for (j=ys; j<ye-2; j++)
            for (i=xs; i<xe-2; i++) {
              double vortx = vortsum[k+1][j+1][i+1].x/N;
              x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( vort2sum[k+1][j+1][i+1].x/N - vortx*vortx );
            }
        if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        I = TECDAT100(&III, x, &DIsDouble);

        for (k=zs; k<ze-2; k++)
          for (j=ys; j<ye-2; j++)
            for (i=xs; i<xe-2; i++) {
              double vorty = vortsum[k+1][j+1][i+1].y/N;
              x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( vort2sum[k+1][j+1][i+1].y/N - vorty*vorty );
            }
        if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        I = TECDAT100(&III, x, &DIsDouble);

        for (k=zs; k<ze-2; k++)
          for (j=ys; j<ye-2; j++)
            for (i=xs; i<xe-2; i++) {
              double vortz = vortsum[k+1][j+1][i+1].z/N;
              x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( vort2sum[k+1][j+1][i+1].z/N - vortz*vortz );
            }
        if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        I = TECDAT100(&III, x, &DIsDouble);

        DAVecRestoreArray(/* user[bi]. */ user->fda, Vort_sum, &vortsum);
        DAVecRestoreArray(/* user[bi]. */ user->fda, Vort_square_sum, &vort2sum);

        VecDestroy(Vort_sum);
        VecDestroy(Vort_square_sum);

      }
    }
    else if(avg==2) {
      PetscViewer viewer;
      Vec Conc_sum;
      PetscReal ***conc_sum;
      Vec K_sum;
      PetscReal ***ksum;
      Vec Level_sum;
      Vec Level_square_sum;
      PetscReal ***level_sum;
      PetscReal ***level2sum;
      char filen[128];

//      PetscPrintf(PETSC_COMM_WORLD, "here!\n");



      DACreateGlobalVector(user->fda, &user->Ucat_sum);
      sprintf(filen, "su0_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoadIntoVector(viewer, (user->Ucat_sum));
      PetscViewerDestroy(viewer);
      DAVecGetArray(/* user[bi]. */ user->fda, /* user[bi]. */ user->Ucat_sum, &usum);

      //PetscPrintf(PETSC_COMM_WORLD, "here!\n");

      if(!sandwave){
        DACreateGlobalVector(user->fda, &user->Ucat_square_sum);
        if(rans) {
          DACreateGlobalVector(user->da, &K_sum);
        }

        if(conv_diff) {
          DACreateGlobalVector(user->da, &Conc_sum);
        }


        if(conv_diff) {
          sprintf(filen, "sconc_%06d_%1d.dat", ti, user->_this);
          PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
          VecLoadIntoVector(viewer, Conc_sum);
          PetscViewerDestroy(viewer);
        }

        if(levelset) {
          DACreateGlobalVector(user->da, &Level_sum);
          DACreateGlobalVector(user->da, &Level_square_sum);
        }

        if(levelset) {
          sprintf(filen, "slevel_%06d_%1d.dat", ti, user->_this);
          PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
          VecLoadIntoVector(viewer, Level_sum);
          PetscViewerDestroy(viewer);

          sprintf(filen, "slevel2_%06d_%1d.dat", ti, user->_this);
          PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
          VecLoadIntoVector(viewer, Level_square_sum);
          PetscViewerDestroy(viewer);
        }

        if(rans) {
          sprintf(filen, "sk_%06d_%1d.dat", ti, user->_this);
          PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
          VecLoadIntoVector(viewer, K_sum);
          PetscViewerDestroy(viewer);
        }
        else {
          sprintf(filen, "su2_%06d_%1d.dat", ti, user->_this);
          PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
          VecLoadIntoVector(viewer, (user->Ucat_square_sum));
          PetscViewerDestroy(viewer);
        }


        if(conv_diff) DAVecGetArray(/* user[bi]. */ user->da, Conc_sum, &conc_sum);

        if(levelset) DAVecGetArray(/* user[bi]. */ user->da, Level_sum, &level_sum);
        if(levelset) DAVecGetArray(/* user[bi]. */ user->da, Level_square_sum, &level2sum);

        if(rans) DAVecGetArray(/* user[bi]. */ user->da, K_sum, &ksum);
        else DAVecGetArray(/* user[bi]. */ user->fda, /* user[bi]. */ user->Ucat_square_sum, &u2sum);

        // U
        for (k=zs; k<ze-2; k++)
          for (j=ys; j<ye-2; j++)
            for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = usum[k+1][j+1][i+1].x/N;
        if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        I = TECDAT100(&III, x, &DIsDouble);

        // V
        for (k=zs; k<ze-2; k++)
          for (j=ys; j<ye-2; j++)
            for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = usum[k+1][j+1][i+1].y/N;
        if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        I = TECDAT100(&III, x, &DIsDouble);

        // W
        for (k=zs; k<ze-2; k++)
          for (j=ys; j<ye-2; j++)
            for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = usum[k+1][j+1][i+1].z/N;
        if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        I = TECDAT100(&III, x, &DIsDouble);

        // k
        for (k=zs; k<ze-2; k++)
          for (j=ys; j<ye-2; j++)
            for(i=xs; i<xe-2; i++) {
              if(rans)  {
                x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ksum[k+1][j+1][i+1]/N;
              }
              else {
                double U = usum[k+1][j+1][i+1].x/N;
                double V = usum[k+1][j+1][i+1].y/N;
                double W = usum[k+1][j+1][i+1].z/N;
                x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ( u2sum[k+1][j+1][i+1].x/N - U*U ) + ( u2sum[k+1][j+1][i+1].y/N - V*V ) + ( u2sum[k+1][j+1][i+1].z/N - W*W );
                x[k * (mx-2)*(my-2) + j*(mx-2) + i] *= 0.5;
              }
            }
        if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        I = TECDAT100(&III, x, &DIsDouble);

        // Conc
        if(conv_diff){
          for (k=zs; k<ze-2; k++)
            for (j=ys; j<ye-2; j++)
              for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = conc_sum[k+1][j+1][i+1]/N;
          if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
          if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
          if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
          if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
          if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
          I = TECDAT100(&III, x, &DIsDouble);
        }

        // Level
        if(levelset){
          for (k=zs; k<ze-2; k++)
            for (j=ys; j<ye-2; j++)
              for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = level_sum[k+1][j+1][i+1]/N;
          if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
          if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
          if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
          if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
          if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
          I = TECDAT100(&III, x, &DIsDouble);

          for (k=zs; k<ze-2; k++)
            for (j=ys; j<ye-2; j++)
              for (i=xs; i<xe-2; i++){
                double L = level_sum[k+1][j+1][i+1]/N;
                x[k * (mx-2)*(my-2) + j*(mx-2) + i] = sqrt( level2sum[k+1][j+1][i+1]/N - L*L ); // rms
              }
          if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
          if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
          if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
          if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
          if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
          I = TECDAT100(&III, x, &DIsDouble);
        }

        if(conv_diff) DAVecRestoreArray(/* user[bi]. */ user->da, Conc_sum, &conc_sum);
        if(levelset) DAVecRestoreArray(/* user[bi]. */ user->da, Level_sum, &level_sum);
        if(levelset) DAVecRestoreArray(/* user[bi]. */ user->da, Level_square_sum, &level2sum);
        if(rans) DAVecRestoreArray(/* user[bi]. */ user->da, K_sum, &ksum);
        else DAVecRestoreArray(/* user[bi]. */ user->fda, /* user[bi]. */ user->Ucat_square_sum, &u2sum);

        if(conv_diff) VecDestroy(Conc_sum);
        if(levelset) VecDestroy(Level_sum);
        if(levelset) VecDestroy(Level_square_sum);
        if(rans) VecDestroy(K_sum);
        else VecDestroy(user->Ucat_square_sum);
      }
      if(sandwave){

        DACreateGlobalVector(/* user[bi]. */ user->fda, &/* user[bi]. */ user->Ucat);
        sprintf(filen, "ufield%06d_%1d.dat", ti, /* user[bi]. */ user->_this);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
        VecLoadIntoVector(viewer, (/* user[bi]. */ user->Ucat));
        PetscViewerDestroy(viewer);
        DAVecGetArray(/* user[bi]. */ user->fda, /* user[bi]. */ user->Ucat, &ucat);

        PetscReal ***ustar_;
        Vec lUstar_;
        //PetscViewer   viewer;
        DACreateGlobalVector(user->da, &lUstar_);
        sprintf(filen, "ustarfield%06d_%1d.dat", ti, user->_this);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
        VecLoadIntoVector(viewer, lUstar_);
        PetscViewerDestroy(viewer);
        DAVecGetArray(/* user[bi]. */ user->da, lUstar_, &ustar_);

        PetscReal ***ustar_sum;
        Vec lUstar_sum;
        //PetscViewer   viewer;
        DACreateGlobalVector(user->da, &lUstar_sum);
        sprintf(filen, "sustar_%06d_%1d.dat", ti, user->_this);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
        VecLoadIntoVector(viewer, lUstar_sum);
        PetscViewerDestroy(viewer);
        DAVecGetArray(/* user[bi]. */ user->da, lUstar_sum, &ustar_sum);



        // u-fluctuations
        for (k=zs; k<ze-2; k++)
          for (j=ys; j<ye-2; j++)
            for (i=xs; i<xe-2; i++) {x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ucat[k+1][j+1][i+1].x-usum[k+1][j+1][i+1].x/N;}
        if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        I = TECDAT100(&III, x, &DIsDouble);

        // v-fluctuations
        for (k=zs; k<ze-2; k++)
          for (j=ys; j<ye-2; j++)
            for (i=xs; i<xe-2; i++) {x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ucat[k+1][j+1][i+1].y-usum[k+1][j+1][i+1].y/N;}
        if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        I = TECDAT100(&III, x, &DIsDouble);

        // w-fluctuations
        for (k=zs; k<ze-2; k++)
          for (j=ys; j<ye-2; j++)
            for (i=xs; i<xe-2; i++) {x[k * (mx-2)*(my-2) + j*(mx-2) + i] = ucat[k+1][j+1][i+1].z-usum[k+1][j+1][i+1].z/N;}
        if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        I = TECDAT100(&III, x, &DIsDouble);

        // F-turbulent sweep
        for (k=zs; k<ze-2; k++)
          for (j=ys; j<ye-2; j++)
            for(i=xs; i<xe-2; i++) {double w_prime = ucat[k+1][j+1][i+1].z-usum[k+1][j+1][i+1].z/N;
              double u_prime = ucat[k+1][j+1][i+1].x-usum[k+1][j+1][i+1].x/N;
              if(u_prime>0.0 && w_prime<0.0) {x[k * (mx-2)*(my-2) + j*(mx-2) + i] = 1.000000;}  else {x[k * (mx-2)*(my-2) + j*(mx-2) + i] = 0.000000;}}
        if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        I = TECDAT100(&III, x, &DIsDouble);


        // ustar
        for (k=zs; k<ze-2; k++)
          for (j=ys; j<ye-2; j++)
            for (i=xs; i<xe-2; i++) {
              x[k * (mx-2)*(my-2) + j*(mx-2) + i] =  ustar_[k+1][j+1][i+1];
            }
        if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        I = TECDAT100(&III, x, &DIsDouble);

        // ustar-fluctuations
        for (k=zs; k<ze-2; k++)
          for (j=ys; j<ye-2; j++)
            for (i=xs; i<xe-2; i++) {
              x[k * (mx-2)*(my-2) + j*(mx-2) + i] =  ustar_[k+1][j+1][i+1]-ustar_sum[k+1][j+1][i+1]/N;
            }
        if(i_average) Iavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(j_average) Javg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(k_average) Kavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ik_average) IKavg(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        if(ikc_average) IKavg_c(x, xs, xe, ys, ye, zs, ze, mx, my, mz);
        I = TECDAT100(&III, x, &DIsDouble);

        DAVecRestoreArray(/* user[bi]. */ user->fda, /* user[bi]. */ user->Ucat, &ucat);
        DAVecRestoreArray(/* user[bi]. */ user->da, lUstar_, &ustar_);
        DAVecRestoreArray(/* user[bi]. */ user->da, lUstar_sum, &ustar_sum);
        VecDestroy(lUstar_);
        VecDestroy(lUstar_sum);
      }
      DAVecRestoreArray(/* user[bi]. */ user->fda, /* user[bi]. */ user->Ucat_sum, &usum);
      VecDestroy(user->Ucat_sum);

    }



    DAVecGetArray(/* user[bi]. */ user->da, /* user[bi]. */ user->Nvert, &nvert);
    for (k=zs; k<ze-2; k++)
      for (j=ys; j<ye-2; j++)
        for (i=xs; i<xe-2; i++) x[k * (mx-2)*(my-2) + j*(mx-2) + i] = nvert[k+1][j+1][i+1];
    I = TECDAT100(&III, x, &DIsDouble);
    DAVecRestoreArray(/* user[bi]. */ user->da, /* user[bi]. */ user->Nvert, &nvert);

    delete []x;
  }
  I = TECEND100();
  return 0;
}

PetscErrorCode TECIOOutQ(std::shared_ptr< UserCtx > & user, int Q)
{
  PetscInt bi;

  char filen[80];
  sprintf(filen, "QCriteria%06d.plt", ti);

  INTEGER4 I, Debug, VIsDouble, DIsDouble, III, IMax, JMax, KMax;
  VIsDouble = 0;
  DIsDouble = 0;
  Debug = 0;

  if(Q==1) {
    printf("qcr=%d, Q-Criterion\n", Q);
    //I = TECINI100((char*)"Result", (char*)"X Y Z Q Velocity_Magnitude Nv",   filen,  (char*)".",   &Debug,  &VIsDouble);
    I = TECINI100((char*)"Result", (char*)"X Y Z Q u v w",   filen,  (char*)".",   &Debug,  &VIsDouble);
    //I = TECINI100((char*)"Result", (char*)"X Y Z Q Velocity_Magnitude",   filen,  (char*)".",   &Debug,  &VIsDouble);
  }
  else if(Q==4) {
    printf("qcr=%d, Q-Criterion\n", Q);
    //I = TECINI100((char*)"Result", (char*)"X Y Z Q Velocity_Magnitude Nv",   filen,  (char*)".",   &Debug,  &VIsDouble);
    I = TECINI100((char*)"Result", (char*)"X Y Z U V W Q",   filen,  (char*)".",   &Debug,  &VIsDouble);
  }
  else if(Q==2) {
    printf("Lambda2-Criterion\n");
    //I = TECINI100((char*)"Result", (char*)"X Y Z Lambda2 Velocity_Magnitude Nv",   filen,  (char*)".",   &Debug,  &VIsDouble);
    I = TECINI100((char*)"Result", (char*)"X Y Z Lambda2 Velocity_Magnitude",   filen,  (char*)".",   &Debug,  &VIsDouble);
  }
  else if(Q==3) {
    printf("Q-Criterion from saved file\n");
    I = TECINI100((char*)"Result", (char*)"X Y Z Q Velocity_Magnitude",   filen,  (char*)".",   &Debug,  &VIsDouble);
  }

  for (bi=0; bi<block_number; bi++) {
    DA da = /* user[bi]. */ user->da, fda = /* user[bi]. */ user->fda;
    DALocalInfo info = /* user[bi]. */ user->info;

    PetscInt        xs = info.xs, xe = info.xs + info.xm;
    PetscInt        ys = info.ys, ye = info.ys + info.ym;
    PetscInt        zs = info.zs, ze = info.zs + info.zm;
    PetscInt        mx = info.mx, my = info.my, mz = info.mz;

    PetscInt        lxs, lys, lzs, lxe, lye, lze;
    PetscInt        i, j, k;
    PetscReal       ***aj;
    Cmpnts  ***ucat, ***coor, ***ucat_o;
    PetscReal       ***p, ***nvert;
    Vec             Coor, zCoor, nCoor;
    VecScatter      ctx;

    Vec X, Y, Z, U, V, W;

    INTEGER4        ZoneType=0, ICellMax=0, JCellMax=0, KCellMax=0;
    INTEGER4        IsBlock=1, NumFaceConnections=0, FaceNeighborMode=0;
    INTEGER4    ShareConnectivityFromZone=0;
    INTEGER4        LOC[8] = {1, 1, 1, 0, 0, 0, 0}; /* 1 is cell-centered 0 is node centered */

    IMax = mx-1; JMax = my-1; KMax = mz-1;

    I = TECZNE100((char*)"Block 1",
                  &ZoneType,      /* Ordered zone */
                  &IMax,
                  &JMax,
                  &KMax,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &IsBlock,       /* ISBLOCK  BLOCK format */
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  LOC,
                  NULL,
                  &ShareConnectivityFromZone); /* No connectivity sharing */

    III = (mx-1) * (my-1) * (mz-1);
    float   *x;
    PetscMalloc(mx*my*mz*sizeof(float), &x);        // seokkoo

    DAGetCoordinates(da, &Coor);
    DAVecGetArray(fda, Coor, &coor);

    for (k=zs; k<ze-1; k++)
      for (j=ys; j<ye-1; j++)
        for (i=xs; i<xe-1; i++) {
          x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].x;
        }
    I = TECDAT100(&III, x, &DIsDouble);

    for (k=zs; k<ze-1; k++)
      for (j=ys; j<ye-1; j++)
        for (i=xs; i<xe-1; i++) {
          x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].y;
        }
    I = TECDAT100(&III, x, &DIsDouble);

    for (k=zs; k<ze-1; k++)
      for (j=ys; j<ye-1; j++)
        for (i=xs; i<xe-1; i++) {
          x[k * (mx-1)*(my-1) + j*(mx-1) + i] = coor[k][j][i].z;
        }

    I = TECDAT100(&III, x, &DIsDouble);
    DAVecRestoreArray(fda, Coor, &coor);

    III = (mx-2) * (my-2) * (mz-2);


    if(Q==1 || Q==4) {
      QCriteria(user);
      DAVecGetArray(/* user[bi]. */ user->da, /* user[bi]. */ user->P, &p);
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) {
            x[k * (mx-2)*(my-2) + j*(mx-2) + i] =p[k+1][j+1][i+1];
          }
      I = TECDAT100(&III, x, &DIsDouble);
      DAVecRestoreArray(/* user[bi]. */ user->da, /* user[bi]. */ user->P, &p);
      //PetscPrintf(PETSC_COMM_WORLD, "here!\n");
    }
    else if(Q==2) {
      Lambda2(user);
      DAVecGetArray(/* user[bi]. */ user->da, /* user[bi]. */ user->P, &p);
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) {
            x[k * (mx-2)*(my-2) + j*(mx-2) + i] = p[k+1][j+1][i+1];
          }
      I = TECDAT100(&III, x, &DIsDouble);
      DAVecRestoreArray(/* user[bi]. */ user->da, /* user[bi]. */ user->P, &p);
    }
    else if(Q==3) {
      char filen2[128];
      PetscViewer     viewer;

      sprintf(filen2, "qfield%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen2, FILE_MODE_READ, &viewer);
      VecLoadIntoVector(viewer, (/* user[bi]. */ user->P));
      PetscViewerDestroy(viewer);

      DAVecGetArray(/* user[bi]. */ user->da, /* user[bi]. */ user->P, &p);
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) {
            x[k * (mx-2)*(my-2) + j*(mx-2) + i] = p[k+1][j+1][i+1];
          }
      I = TECDAT100(&III, x, &DIsDouble);
      DAVecRestoreArray(/* user[bi]. */ user->da, /* user[bi]. */ user->P, &p);
    }
    //if(Q==1 || Q==2 || Q==3){
    if( Q==2 || Q==3){
      Velocity_Magnitude(user);
      DAVecGetArray(/* user[bi]. */ user->da, /* user[bi]. */ user->P, &p);
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) {
            x[k * (mx-2)*(my-2) + j*(mx-2) + i] = p[k+1][j+1][i+1];
          }
      I = TECDAT100(&III, x, &DIsDouble);
      DAVecRestoreArray(/* user[bi]. */ user->da, /* user[bi]. */ user->P, &p);
    }

    //PetscPrintf(PETSC_COMM_WORLD, "here!\n");
    PetscInt iiii;
    if(Q==1){
      iiii=1;
      Velocity_Magnitude_4(user, iiii);
      DAVecGetArray(/* user[bi]. */ user->da, /* user[bi]. */ user->P, &p);
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) {
            x[k * (mx-2)*(my-2) + j*(mx-2) + i] = p[k+1][j+1][i+1];
          }
      I = TECDAT100(&III, x, &DIsDouble);
      DAVecRestoreArray(/* user[bi]. */ user->da, /* user[bi]. */ user->P, &p);
    }
//PetscPrintf(PETSC_COMM_WORLD, "here!\n");
    if(Q==1){
      iiii=2;
      Velocity_Magnitude_4(user, iiii);
      DAVecGetArray(/* user[bi]. */ user->da, /* user[bi]. */ user->P, &p);
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) {
            x[k * (mx-2)*(my-2) + j*(mx-2) + i] = p[k+1][j+1][i+1];
          }
      I = TECDAT100(&III, x, &DIsDouble);
      DAVecRestoreArray(/* user[bi]. */ user->da, /* user[bi]. */ user->P, &p);
    }
    //PetscPrintf(PETSC_COMM_WORLD, "here!\n");
    if(Q==1){
      iiii=3;
      Velocity_Magnitude_4(user, iiii);
      DAVecGetArray(/* user[bi]. */ user->da, /* user[bi]. */ user->P, &p);
      for (k=zs; k<ze-2; k++)
        for (j=ys; j<ye-2; j++)
          for (i=xs; i<xe-2; i++) {
            x[k * (mx-2)*(my-2) + j*(mx-2) + i] = p[k+1][j+1][i+1];
          }
      I = TECDAT100(&III, x, &DIsDouble);
      DAVecRestoreArray(/* user[bi]. */ user->da, /* user[bi]. */ user->P, &p);
    }


    // DAVecGetArray(/* user[bi]. */ user->da, /* user[bi]. */ user->Nvert, &nvert);
    // for (k=zs; k<ze-2; k++)
    // for (j=ys; j<ye-2; j++)
    // for (i=xs; i<xe-2; i++) {
    // x[k * (mx-2)*(my-2) + j*(mx-2) + i] = nvert[k+1][j+1][i+1];
    // }
    // I = TECDAT100(&III, x, &DIsDouble);
    // DAVecRestoreArray(/* user[bi]. */ user->da, /* user[bi]. */ user->Nvert, &nvert);

    PetscFree(x);
  }
  I = TECEND100();

  //PetscPrintf(PETSC_COMM_WORLD, "last here!\n");
  return 0;
}

PetscErrorCode FormMetrics(std::shared_ptr< UserCtx > & user)
{
  DA            cda;
  Cmpnts        ***csi, ***eta, ***zet;
  PetscScalar   ***aj;
  Vec           coords;
  Cmpnts        ***coor;

  DA            da = user->da, fda = user->fda;
  Vec           Csi = user->Csi, Eta = user->Eta, Zet = user->Zet;
  Vec           Aj = user->Aj;
  Vec           ICsi = user->ICsi, IEta = user->IEta, IZet = user->IZet;
  Vec           JCsi = user->JCsi, JEta = user->JEta, JZet = user->JZet;
  Vec           KCsi = user->KCsi, KEta = user->KEta, KZet = user->KZet;
  Vec           IAj = user->IAj, JAj = user->JAj, KAj = user->KAj;


  Cmpnts        ***icsi, ***ieta, ***izet;
  Cmpnts        ***jcsi, ***jeta, ***jzet;
  Cmpnts        ***kcsi, ***keta, ***kzet;
  Cmpnts        ***gs;
  PetscReal     ***iaj, ***jaj, ***kaj;

  Vec           Cent = user->Cent; //local working array for storing cell center geometry

  Vec           Centx, Centy, Centz, lCoor;
  Cmpnts        ***cent, ***centx, ***centy, ***centz;

  PetscInt      xs, ys, zs, xe, ye, ze;
  DALocalInfo   info;

  PetscInt      mx, my, mz;
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  PetscScalar   dxdc, dydc, dzdc, dxde, dyde, dzde, dxdz, dydz, dzdz;
  PetscInt      i, j, k, ia, ja, ka, ib, jb, kb;
  PetscInt      gxs, gxe, gys, gye, gzs, gze;
  PetscErrorCode        ierr;

  PetscReal     xcp, ycp, zcp, xcm, ycm, zcm;
  DAGetLocalInfo(da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  DAGetCoordinateDA(da, &cda);
  DAVecGetArray(cda, Csi, &csi);
  DAVecGetArray(cda, Eta, &eta);
  DAVecGetArray(cda, Zet, &zet);
  ierr = DAVecGetArray(da, Aj,  &aj); CHKERRQ(ierr);

  DAGetGhostedCoordinates(da, &coords);
  DAVecGetArray(fda, coords, &coor);


  //  VecDuplicate(coords, &Cent);
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  /* Calculating transformation metrics in i direction */
  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=xs; i<lxe; i++){
        /* csi = de X dz */
        dxde = 0.5 * (coor[k  ][j  ][i  ].x + coor[k-1][j  ][i  ].x -
                      coor[k  ][j-1][i  ].x - coor[k-1][j-1][i  ].x);
        dyde = 0.5 * (coor[k  ][j  ][i  ].y + coor[k-1][j  ][i  ].y -
                      coor[k  ][j-1][i  ].y - coor[k-1][j-1][i  ].y);
        dzde = 0.5 * (coor[k  ][j  ][i  ].z + coor[k-1][j  ][i  ].z -
                      coor[k  ][j-1][i  ].z - coor[k-1][j-1][i  ].z);

        dxdz = 0.5 * (coor[k  ][j-1][i  ].x + coor[k  ][j  ][i  ].x -
                      coor[k-1][j-1][i  ].x - coor[k-1][j  ][i  ].x);
        dydz = 0.5 * (coor[k  ][j-1][i  ].y + coor[k  ][j  ][i  ].y -
                      coor[k-1][j-1][i  ].y - coor[k-1][j  ][i  ].y);
        dzdz = 0.5 * (coor[k  ][j-1][i  ].z + coor[k  ][j  ][i  ].z -
                      coor[k-1][j-1][i  ].z - coor[k-1][j  ][i  ].z);

        csi[k][j][i].x = dyde * dzdz - dzde * dydz;
        csi[k][j][i].y =-dxde * dzdz + dzde * dxdz;
        csi[k][j][i].z = dxde * dydz - dyde * dxdz;


      }
    }
  }

  // Need more work -- lg65
  /* calculating j direction metrics */
  for (k=lzs; k<lze; k++){
    for (j=ys; j<lye; j++){
      for (i=lxs; i<lxe; i++){

        /* eta = dz X de */
        dxdc = 0.5 * (coor[k  ][j  ][i  ].x + coor[k-1][j  ][i  ].x -
                      coor[k  ][j  ][i-1].x - coor[k-1][j  ][i-1].x);
        dydc = 0.5 * (coor[k  ][j  ][i  ].y + coor[k-1][j  ][i  ].y -
                      coor[k  ][j  ][i-1].y - coor[k-1][j  ][i-1].y);
        dzdc = 0.5 * (coor[k  ][j  ][i  ].z + coor[k-1][j  ][i  ].z -
                      coor[k  ][j  ][i-1].z - coor[k-1][j  ][i-1].z);

        dxdz = 0.5 * (coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x -
                      coor[k-1][j  ][i  ].x - coor[k-1][j  ][i-1].x);
        dydz = 0.5 * (coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y -
                      coor[k-1][j  ][i  ].y - coor[k-1][j  ][i-1].y);
        dzdz = 0.5 * (coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z -
                      coor[k-1][j  ][i  ].z - coor[k-1][j  ][i-1].z);

        eta[k][j][i].x = dydz * dzdc - dzdz * dydc;
        eta[k][j][i].y =-dxdz * dzdc + dzdz * dxdc;
        eta[k][j][i].z = dxdz * dydc - dydz * dxdc;

      }
    }
  }

  /* calculating k direction metrics */
  for (k=zs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){
        dxdc = 0.5 * (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x -
                      coor[k  ][j  ][i-1].x - coor[k  ][j-1][i-1].x);
        dydc = 0.5 * (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y -
                      coor[k  ][j  ][i-1].y - coor[k  ][j-1][i-1].y);
        dzdc = 0.5 * (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z -
                      coor[k  ][j  ][i-1].z - coor[k  ][j-1][i-1].z);

        dxde = 0.5 * (coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x -
                      coor[k  ][j-1][i  ].x - coor[k  ][j-1][i-1].x);
        dyde = 0.5 * (coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y -
                      coor[k  ][j-1][i  ].y - coor[k  ][j-1][i-1].y);
        dzde = 0.5 * (coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z -
                      coor[k  ][j-1][i  ].z - coor[k  ][j-1][i-1].z);

        zet[k][j][i].x = dydc * dzde - dzdc * dyde;
        zet[k][j][i].y =-dxdc * dzde + dzdc * dxde;
        zet[k][j][i].z = dxdc * dyde - dydc * dxde;

        /*      if ((i==1 || i==mx-2) && j==1 && (k==1 || k==0)) {
                PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", dxdc * dyde, dydc * dxde, dzdc);
                PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", dxde, dyde, dzde);
                PetscPrintf(PETSC_COMM_WORLD, "Met %e %e %e\n", zet[k][j][i].x, zet[k][j][i].y, zet[k][j][i].z);
                }*/

      }
    }
  }

  /* calculating Jacobian of the transformation */
  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){

        if (i>0 && j>0 && k>0) {
          dxdc = 0.25 * (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
                         coor[k-1][j  ][i  ].x + coor[k-1][j-1][i  ].x -
                         coor[k  ][j  ][i-1].x - coor[k  ][j-1][i-1].x -
                         coor[k-1][j  ][i-1].x - coor[k-1][j-1][i-1].x);
          dydc = 0.25 * (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
                         coor[k-1][j  ][i  ].y + coor[k-1][j-1][i  ].y -
                         coor[k  ][j  ][i-1].y - coor[k  ][j-1][i-1].y -
                         coor[k-1][j  ][i-1].y - coor[k-1][j-1][i-1].y);
          dzdc = 0.25 * (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
                         coor[k-1][j  ][i  ].z + coor[k-1][j-1][i  ].z -
                         coor[k  ][j  ][i-1].z - coor[k  ][j-1][i-1].z -
                         coor[k-1][j  ][i-1].z - coor[k-1][j-1][i-1].z);

          dxde = 0.25 * (coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x +
                         coor[k-1][j  ][i  ].x + coor[k-1][j  ][i-1].x -
                         coor[k  ][j-1][i  ].x - coor[k  ][j-1][i-1].x -
                         coor[k-1][j-1][i  ].x - coor[k-1][j-1][i-1].x);
          dyde = 0.25 * (coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y +
                         coor[k-1][j  ][i  ].y + coor[k-1][j  ][i-1].y -
                         coor[k  ][j-1][i  ].y - coor[k  ][j-1][i-1].y -
                         coor[k-1][j-1][i  ].y - coor[k-1][j-1][i-1].y);
          dzde = 0.25 * (coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z +
                         coor[k-1][j  ][i  ].z + coor[k-1][j  ][i-1].z -
                         coor[k  ][j-1][i  ].z - coor[k  ][j-1][i-1].z -
                         coor[k-1][j-1][i  ].z - coor[k-1][j-1][i-1].z);

          dxdz = 0.25 * (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
                         coor[k  ][j  ][i-1].x + coor[k  ][j-1][i-1].x -
                         coor[k-1][j  ][i  ].x - coor[k-1][j-1][i  ].x -
                         coor[k-1][j  ][i-1].x - coor[k-1][j-1][i-1].x);
          dydz = 0.25 * (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
                         coor[k  ][j  ][i-1].y + coor[k  ][j-1][i-1].y -
                         coor[k-1][j  ][i  ].y - coor[k-1][j-1][i  ].y -
                         coor[k-1][j  ][i-1].y - coor[k-1][j-1][i-1].y);
          dzdz = 0.25 * (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
                         coor[k  ][j  ][i-1].z + coor[k  ][j-1][i-1].z -
                         coor[k-1][j  ][i  ].z - coor[k-1][j-1][i  ].z -
                         coor[k-1][j  ][i-1].z - coor[k-1][j-1][i-1].z);

          aj[k][j][i] = dxdc * (dyde * dzdz - dzde * dydz) -
            dydc * (dxde * dzdz - dzde * dxdz) +
            dzdc * (dxde * dydz - dyde * dxdz);
          aj[k][j][i] = 1./aj[k][j][i];

#ifdef NEWMETRIC
          csi[k][j][i].x = dyde * dzdz - dzde * dydz;
          csi[k][j][i].y =-dxde * dzdz + dzde * dxdz;
          csi[k][j][i].z = dxde * dydz - dyde * dxdz;

          eta[k][j][i].x = dydz * dzdc - dzdz * dydc;
          eta[k][j][i].y =-dxdz * dzdc + dzdz * dxdc;
          eta[k][j][i].z = dxdz * dydc - dydz * dxdc;

          zet[k][j][i].x = dydc * dzde - dzdc * dyde;
          zet[k][j][i].y =-dxdc * dzde + dzdc * dxde;
          zet[k][j][i].z = dxdc * dyde - dydc * dxde;
#endif
        }
      }
    }
  }

  // mirror grid outside the boundary
  if (xs==0) {
    i = xs;
    for (k=zs; k<ze; k++)
      for (j=ys; j<ye; j++) {
#ifdef NEWMETRIC
        csi[k][j][i] = csi[k][j][i+1];
#endif
        eta[k][j][i] = eta[k][j][i+1];
        zet[k][j][i] = zet[k][j][i+1];
        aj[k][j][i] = aj[k][j][i+1];
      }
  }

  if (xe==mx) {
    i = xe-1;
    for (k=zs; k<ze; k++)
      for (j=ys; j<ye; j++) {
#ifdef NEWMETRIC
        csi[k][j][i] = csi[k][j][i-1];
#endif
        eta[k][j][i] = eta[k][j][i-1];
        zet[k][j][i] = zet[k][j][i-1];
        aj[k][j][i] = aj[k][j][i-1];
      }
  }


  if (ys==0) {
    j = ys;
    for (k=zs; k<ze; k++)
      for (i=xs; i<xe; i++) {
#ifdef NEWMETRIC
        eta[k][j][i] = eta[k][j+1][i];
#endif
        csi[k][j][i] = csi[k][j+1][i];
        zet[k][j][i] = zet[k][j+1][i];
        aj[k][j][i] = aj[k][j+1][i];
      }
  }


  if (ye==my) {
    j = ye-1;
    for (k=zs; k<ze; k++)
      for (i=xs; i<xe; i++) {
#ifdef NEWMETRIC
        eta[k][j][i] = eta[k][j-1][i];
#endif
        csi[k][j][i] = csi[k][j-1][i];
        zet[k][j][i] = zet[k][j-1][i];
        aj[k][j][i] = aj[k][j-1][i];
      }
  }


  if (zs==0) {
    k = zs;
    for (j=ys; j<ye; j++)
      for (i=xs; i<xe; i++) {
#ifdef NEWMETRIC
        zet[k][j][i] = zet[k+1][j][i];
#endif
        eta[k][j][i] = eta[k+1][j][i];
        csi[k][j][i] = csi[k+1][j][i];
        aj[k][j][i] = aj[k+1][j][i];
      }
  }


  if (ze==mz) {
    k = ze-1;
    for (j=ys; j<ye; j++)
      for (i=xs; i<xe; i++) {
#ifdef NEWMETRIC
        zet[k][j][i] = zet[k-1][j][i];
#endif
        eta[k][j][i] = eta[k-1][j][i];
        csi[k][j][i] = csi[k-1][j][i];
        aj[k][j][i] = aj[k-1][j][i];
      }
  }


  //  PetscPrintf(PETSC_COMM_WORLD, "Local info: %d", info.mx);



  DAVecRestoreArray(cda, Csi, &csi);
  DAVecRestoreArray(cda, Eta, &eta);
  DAVecRestoreArray(cda, Zet, &zet);
  DAVecRestoreArray(da, Aj,  &aj);


  DAVecRestoreArray(cda, coords, &coor);


  VecAssemblyBegin(Csi);
  VecAssemblyEnd(Csi);
  VecAssemblyBegin(Eta);
  VecAssemblyEnd(Eta);
  VecAssemblyBegin(Zet);
  VecAssemblyEnd(Zet);
  VecAssemblyBegin(Aj);
  VecAssemblyEnd(Aj);

  PetscBarrier(PETSC_NULL);
  return 0;
}

PetscErrorCode Ucont_P_Binary_Input(std::shared_ptr< UserCtx > & user)
{
  PetscViewer   viewer;

  char filen2[128];

  PetscOptionsClearValue("-vecload_block_size");
  sprintf(filen2, "pfield%06d_%1d.dat", ti, user->_this);

  PetscViewer   pviewer;
  //Vec temp;
  PetscInt rank;
  PetscReal norm;

  if(file_exist(filen2))
    if(!onlyV) {
      //DACreateNaturalVector(user->da, &temp);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen2, FILE_MODE_READ, &pviewer);
      VecLoadIntoVector(pviewer, (user->P));
      VecNorm(user->P, NORM_INFINITY, &norm);
      PetscPrintf(PETSC_COMM_WORLD, "PIn %le\n", norm);
      PetscViewerDestroy(pviewer);
      //VecDestroy(temp);
    }

/*  if(conv_diff){
    sprintf(filen2, "cfield%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen2, FILE_MODE_READ, &pviewer);
    VecLoadIntoVector(pviewer, (user->Conc));
    PetscViewerDestroy(pviewer);
    }
*/
  if(nv_once) sprintf(filen2, "nvfield%06d_%1d.dat", 0, user->_this);
  else sprintf(filen2, "nvfield%06d_%1d.dat", ti, user->_this);

  if(cs) sprintf(filen2, "cs_%06d_%1d.dat", ti, user->_this);

  if( !nv_once || (nv_once && ti==tis) )
  {
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen2, FILE_MODE_READ, &pviewer);
    VecLoadIntoVector(pviewer, (user->Nvert));
    PetscViewerDestroy(pviewer);
  }

  return false;
}

PetscErrorCode Ucont_P_Binary_Input1(std::shared_ptr< UserCtx > & user)
{
  PetscViewer viewer;
  char filen[128];

  sprintf(filen, "ufield%06d_%1d.dat", ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);

  PetscInt N;

  VecGetSize(user->Ucat, &N);
  PetscPrintf(PETSC_COMM_WORLD, "PPP %d\n", N);
  VecLoadIntoVector(viewer, (user->Ucat));
  PetscViewerDestroy(viewer);

  PetscBarrier(PETSC_NULL);

  return false;
}

PetscErrorCode Ucont_P_Binary_Input_Averaging(std::shared_ptr< UserCtx > & user)
{
  PetscViewer viewer;
  char filen[128];
  /*
    sprintf(filen, "su0_%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoadIntoVector(viewer, (user->Ucat_sum));
    PetscViewerDestroy(viewer);

    sprintf(filen, "su1_%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoadIntoVector(viewer, (user->Ucat_cross_sum));
    PetscViewerDestroy(viewer);

    sprintf(filen, "su2_%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoadIntoVector(viewer, (user->Ucat_square_sum));
    PetscViewerDestroy(viewer);
  */
  /*
    sprintf(filen, "sp_%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoadIntoVector(viewer, (user->P));
    PetscViewerDestroy(viewer);
  */

  if(pcr) {
    Vec Ptmp;
    VecDuplicate(user->P, &Ptmp);

    sprintf(filen, "pfield%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoadIntoVector(viewer, user->P);
    PetscViewerDestroy(viewer);

    sprintf(filen, "sp_%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoadIntoVector(viewer, Ptmp);
    PetscViewerDestroy(viewer);


    VecScale(Ptmp, -1./((double)tis+1.0));
    VecAXPY(user->P, 1., Ptmp);

    VecDestroy(Ptmp);
  }

  if(nv_once) sprintf(filen, "nvfield%06d_%1d.dat", 0, user->_this);
  else sprintf(filen, "nvfield%06d_%1d.dat", ti, user->_this);

  //if(cs) sprintf(filen2, "cs_%06d_%1d.dat", ti, user->_this);

  if( !nv_once || (nv_once && ti==tis) )
  {
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoadIntoVector(viewer, (user->Nvert));
    PetscViewerDestroy(viewer);
  }
  /*
    if( !nv_once || (nv_once && ti==tis) ) {
    sprintf(filen, "nvfield%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
    VecLoadIntoVector(viewer, (user->Nvert));
    PetscViewerDestroy(viewer);
    }
  */
  PetscBarrier(PETSC_NULL);

  return false;
}

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char **argv)
{
  PetscTruth flag;

  DA      da, fda;
  Vec     qn, qnm;
  Vec     c;
  std::shared_ptr< UserCtx > user(new UserCtx);

  PetscErrorCode ierr;

  IBMNodes        ibm, *ibm0, *ibm1;

  PetscInitialize(&argc, &argv, (char *)0, help);

#if 0
  {
    FILE *fd;
    fd = fopen("control.dat", "r");
    if (!fd)
    {
      printf("\nCould not open file control.dat.\n\n");
      exit(1);
    }
    fclose(fd);
  }

  PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);
  
  char tmp_str[256];
  PetscOptionsGetString(PETSC_NULL, "-prefix", tmp_str, 256, &flag);
  if(flag)sprintf(prefix, "%s_", tmp_str);
  else sprintf(prefix, "");

  PetscOptionsGetReal(PETSC_NULL, "-gx", &gravity_x, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-gy", &gravity_y, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-gz", &gravity_z, PETSC_NULL);

  PetscOptionsGetInt(PETSC_NULL, "-vc", &vc, PETSC_NULL);
#endif
  PetscOptionsGetInt(PETSC_NULL, "-binary", &binary_input, &flag);
#if 0
  PetscOptionsGetInt(PETSC_NULL, "-xyz", &xyz_input, &flag);
  PetscOptionsGetInt(PETSC_NULL, "-rans", &rans, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-conv_diff", &conv_diff, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-sandwave", &sandwave, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-ransout", &rans_output, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-levelset", &levelset, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-avg", &avg, &flag);
  PetscOptionsGetInt(PETSC_NULL, "-shear", &shear, &flag);
  PetscOptionsGetInt(PETSC_NULL, "-averaging", &averaging_option, &flag); // from control.dat

  PetscOptionsGetInt(PETSC_NULL, "-cs", &cs, &flag);
  PetscOptionsGetInt(PETSC_NULL, "-i_periodic", &i_periodic, &flag);
  PetscOptionsGetInt(PETSC_NULL, "-j_periodic", &j_periodic, &flag);
  PetscOptionsGetInt(PETSC_NULL, "-k_periodic", &k_periodic, &flag);

  PetscOptionsGetInt(PETSC_NULL, "-ii_periodic", &i_periodic, &flag);
  PetscOptionsGetInt(PETSC_NULL, "-jj_periodic", &j_periodic, &flag);
  PetscOptionsGetInt(PETSC_NULL, "-kk_periodic", &k_periodic, &flag);

  PetscOptionsGetInt(PETSC_NULL, "-nv", &nv_once, &flag);
  printf("nv_once=%d\n", nv_once);

  int QCR = 0;
  PetscOptionsGetInt(PETSC_NULL, "-qcr", &QCR, PETSC_NULL);
#endif

  PetscOptionsGetInt(PETSC_NULL, "-tis", &tis, &flag);
  if (!flag) PetscPrintf(PETSC_COMM_WORLD, "Need the starting number!\n");

  PetscOptionsGetInt(PETSC_NULL, "-tie", &tie, &flag);
  if (!flag) tie = tis;

  PetscOptionsGetInt(PETSC_NULL, "-ts", &tsteps, &flag);
  if (!flag) tsteps = 5; /* Default increasement is 5 */

  // Option to generate data of all points on or out of the outlet BC.
  PetscOptionsGetInt(PETSC_NULL, "-end", &user->end, &flag);
  if (!flag) user->end = 0;

  // Option to show count of points in cells that went through a cross plane
  // perpendicular to the x axis.
  user->xflag = true;
  PetscOptionsGetReal(PETSC_NULL, "-x", &user->x, &flag);
  if (!flag)
    user->xflag = false;

  // Option to show all the points that went through a cross plane
  // perpendicular to the x axis.
  user->xpartflag = true;
  PetscOptionsGetReal(PETSC_NULL, "-xpart", &user->xpart, &flag);
  if (!flag)
    user->xpartflag = false;

  // Option to show all the points projected on a plane perpendicular to
  // the y axis for one time step.
  user->yconc = true;
  PetscOptionsGetReal(PETSC_NULL, "-yconc", &user->y, &flag);
  if (!flag)
    user->yconc = false;

  // Option to show all the points projected on a plane perpendicular to
  // the y axis for all specified time steps.
  user->yavrg = true;
  PetscOptionsGetReal(PETSC_NULL, "-yavrg", &user->y, &flag);
  if (!flag)
    user->yavrg = false;

  // Option to show all the points projected on a plane perpendicular to
  // the z axis for one time step.
  user->zconc = true;
  PetscOptionsGetReal(PETSC_NULL, "-zconc", &user->z, &flag);
  if (!flag)
    user->zconc = false;

  // Option to show all the points projected on a plane perpendicular to
  // the z axis for all specified time steps.
  user->zavrg = true;
  PetscOptionsGetReal(PETSC_NULL, "-zavrg", &user->z, &flag);
  if (!flag)
    user->zavrg = false;

  // Scale for avrg options. NOTE: program also scales by dividing by the number
  // file. Scale is in addition to this.
  user->scaleFlag = true;
  PetscOptionsGetReal(PETSC_NULL, "-scale", &user->scale, &flag);
  if (!flag)
    user->scaleFlag = false;

  // For zavrg, only include particles that are close to the ground to get
  // the concentration on the ground.
  user->groundFlag = true;
  PetscOptionsGetReal(PETSC_NULL, "-ground", &user->ground, &flag);
  if (!flag)
    user->groundFlag = false;

  // Always output stdout when a newline is found.
  int nobuffer = 0;
  PetscOptionsGetInt(PETSC_NULL, "-nobuffer", &nobuffer, &flag);
  if (nobuffer)
    setbuf(stdout, NULL);

  // Output a particle file suitable for tec360, otherwise it will output
  // a tec360 results files.
  PetscOptionsGetInt(PETSC_NULL, "-partfile", &user->partfile, &flag);
  if (!flag)
    user->partfile = 0;

  // Generate results file with only nvfield and ufield files.
  // As include FTLE results it ftle.output file which contains
  //   the ftle filenames to include in results.
  PetscOptionsGetInt(PETSC_NULL, "-uonly", &user->uonly, &flag);
  if (!flag)
    user->uonly = 0;


  
#if 0
  PetscOptionsGetInt(PETSC_NULL, "-onlyV", &onlyV, &flag);
  PetscOptionsGetInt(PETSC_NULL, "-iavg", &i_average, &flag);
  PetscOptionsGetInt(PETSC_NULL, "-javg", &j_average, &flag);
  PetscOptionsGetInt(PETSC_NULL, "-kavg", &k_average, &flag);
  PetscOptionsGetInt(PETSC_NULL, "-ikavg", &ik_average, &flag);
  PetscOptionsGetInt(PETSC_NULL, "-pcr", &pcr, &flag);
  PetscOptionsGetInt(PETSC_NULL, "-reynolds", &reynolds, &flag);
  PetscOptionsGetInt(PETSC_NULL, "-width_ave_instant", &width_ave_instant, PETSC_NULL);

  PetscOptionsGetInt(PETSC_NULL, "-ikcavg", &ikc_average, &flag);
  if(flag) {
    PetscTruth flag1, flag2;
    PetscOptionsGetInt(PETSC_NULL, "-pi", &pi, &flag1);
    PetscOptionsGetInt(PETSC_NULL, "-pk", &pk, &flag2);

    if(!flag1 || !flag2) {
      printf("To use -ikcavg you must set -pi and -pk, which are number of points in i- and k- directions.\n");
      exit(1);
    }
  }

  PetscOptionsGetInt(PETSC_NULL, "-paraview", &paraview, PETSC_NULL);


  if(pcr) avg=1;
  if(i_average) avg=1;
  if(j_average) avg=1;
  if(k_average) avg=1;
  if(ik_average) avg=1;
  if(ikc_average) avg=1;


  if(i_average + j_average + k_average >1) PetscPrintf(PETSC_COMM_WORLD, "Iavg and Javg cannot be set together !! !\n"), exit(1);
#endif

  PetscInt rank, bi;


  PetscMalloc(sizeof(IBMNodes), &ibm0);
  PetscMalloc(sizeof(IBMNodes), &ibm1);

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(xyz_input) {block_number=1;}
  else {
    FILE *fd;
    fd = fopen("grid.dat", "r");
    if (!fd)
    {
      printf("\nCould not open file grid.dat.\n\n");
      exit(1);
    }
    if(binary_input) fread(&block_number, sizeof(int), 1, fd);
    else fscanf(fd, "%i\n", &block_number);

    // This program currently only handles one block. WRO 20211001
    block_number = 1;

    MPI_Bcast(&block_number, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    fclose(fd);
  }

  PetscOptionsGetReal(PETSC_NULL, "-ren", &user->ren, PETSC_NULL);

  ReadCoordinates(user);

#if 0
  PetscPrintf(PETSC_COMM_WORLD, "read coord!\n");

  for (bi=0; bi<block_number; bi++) {
    DACreateGlobalVector(/* user[bi]. */ user->da, &/* user[bi]. */ user->Nvert);
    if(shear) {
      DACreateGlobalVector(/* user[bi]. */ user->fda, &/* user[bi]. */ user->Csi);
      DACreateGlobalVector(/* user[bi]. */ user->fda, &/* user[bi]. */ user->Eta);
      DACreateGlobalVector(/* user[bi]. */ user->fda, &/* user[bi]. */ user->Zet);
      DACreateGlobalVector(/* user[bi]. */ user->da, &/* user[bi]. */ user->Aj);
      FormMetrics(&(user[bi]));

      Calc_avg_shear_stress(&(user[bi]));

      VecDestroy(/* user[bi]. */ user->Csi);
      VecDestroy(/* user[bi]. */ user->Eta);
      VecDestroy(/* user[bi]. */ user->Zet);
      VecDestroy(/* user[bi]. */ user->Aj);
      exit(0);
    }
    else if(!avg) {
      DACreateGlobalVector(/* user[bi]. */ user->da, &/* user[bi]. */ user->P);
      DACreateGlobalVector(/* user[bi]. */ user->fda, &/* user[bi]. */ user->Ucat);
      //if(conv_diff)DACreateGlobalVector(/* user[bi]. */ user->da, &/* user[bi]. */ user->Conc);
      if(!vc) DACreateGlobalVector(/* user[bi]. */ user->fda, &/* user[bi]. */ user->Ucat_o);

      if(QCR) {
        if(QCR==1 || QCR==2) {
          DACreateGlobalVector(/* user[bi]. */ user->fda, &/* user[bi]. */ user->Csi);
          DACreateGlobalVector(/* user[bi]. */ user->fda, &/* user[bi]. */ user->Eta);
          DACreateGlobalVector(/* user[bi]. */ user->fda, &/* user[bi]. */ user->Zet);
          DACreateGlobalVector(/* user[bi]. */ user->da, &/* user[bi]. */ user->Aj);
          FormMetrics(&(user[bi]));
        }
      }
    }
    else {
      if(pcr) {
        DACreateGlobalVector(/* user[bi]. */ user->da, &/* user[bi]. */ user->P);
        DACreateGlobalVector(/* user[bi]. */ user->fda, &/* user[bi]. */ user->Ucat);
      }
      else if(avg==1) {
        /*
          DACreateGlobalVector(/* user[bi]. */ user->fda, &/* user[bi]. */ user->Ucat_sum);
          DACreateGlobalVector(/* user[bi]. */ user->fda, &/* user[bi]. */ user->Ucat_cross_sum);
          DACreateGlobalVector(/* user[bi]. */ user->fda, &/* user[bi]. */ user->Ucat_square_sum);
        */
      }
      else if(avg==2) {       // just compute k
        /*
          DACreateGlobalVector(/* user[bi]. */ user->fda, &/* user[bi]. */ user->Ucat_sum);
          DACreateGlobalVector(/* user[bi]. */ user->fda, &/* user[bi]. */ user->Ucat_square_sum);
        */
      }
    }

  }



  if(avg) {
    if(i_average) PetscPrintf(PETSC_COMM_WORLD, "Averaging in I direction!\n");
    else if(j_average) PetscPrintf(PETSC_COMM_WORLD, "Averaging in J direction!\n");
    else if(k_average) PetscPrintf(PETSC_COMM_WORLD, "Averaging in K direction!\n");
    else if(ik_average) PetscPrintf(PETSC_COMM_WORLD, "Averaging in IK direction!\n");
    else PetscPrintf(PETSC_COMM_WORLD, "Averaging !\n");
    /*
      DACreateGlobalVector(/* user[bi]. */ user->fda, &user->Ucat_sum);
      DACreateGlobalVector(/* user[bi]. */ user->fda, &user->Ucat_cross_sum);
      DACreateGlobalVector(/* user[bi]. */ user->fda, &user->Ucat_square_sum);
    */

  }
#endif

  // WRO
  Vec Coor;  
  DAGetGhostedCoordinates(user->da, &Coor);
  DAVecGetArray(user->fda, Coor, &user->coor);

  Vec Count;
  DACreateGlobalVector(user->fda, &Count);
  DAVecGetArray(user->fda, Count, &user->count);

  new(&user->spatialIndexGrid)
        std::shared_ptr<SpatialIndexGrid>(nullptr, [](SpatialIndexGrid *) { });

  if (user->xflag)
    new(&user->partMap) std::unordered_map < int, CrossX >();

  LOCAL::genRtree(user);

  // Total number of files read for totals options.
  int fileCnt = 0;

  for (ti=tis; ti<=tie; ti+=tsteps)
  {
    if (ti == tis || (!user->yavrg && !user->zavrg))
      clearCount(user);

#if 0
    for (bi=0; bi<block_number; bi++) {
      if(avg) Ucont_P_Binary_Input_Averaging(&user[bi]);
      else {
        Ucont_P_Binary_Input(&user[bi]);
        Ucont_P_Binary_Input1(&user[bi]);
      }
    }

    if (!QCR) {
      if(avg) TECIOOut_Averaging(user);
      else TECIOOut_V(user, onlyV);
      //TECIOOut(user);
    }
    else {
      TECIOOutQ(user, QCR);
    }
#endif

    if (user->xflag || user->xpartflag)
    {
      char filename[80];
      sprintf(filename, "Particle%06d_0.dat", ti);
      printf("Reading file %s\n", filename);
      calcXCrossing(user, ti);
    }
    else if (user->yconc || user->zconc)
    {
      calcCount(user, ti);
      concOutput(user, ti, fileCnt);
    }
    else if (user->yavrg || user->zavrg)
    {
      ++fileCnt;
      calcCount(user, ti);
    }
    else if (user->zmaxflag)
    {
      char filename[80];
      sprintf(filename, "Particle%06d_0.dat", ti);
      printf("\nReading file %s\n", filename);
      calcZMaxCount(user, ti);
    }
    else if (user->uonly)
      outputUandFTLE(user, ti);
    else
    {
      calcPartCount(user, ti);

      INTEGER4 I;
      char filename[80];
      sprintf(filename, "Result%06d.plt", ti);
      printf("\nGenerating file %s\n", filename);

      TECInitialize(user, "X Y Z CNT", filename);
      TECOutputGrid(user);
      TECOutputVector(user, user->count);
      I = TECEND100();
      if (I != 0)
      {
        printf("TECEND100 returned error %d\n", I);
        exit(1);
      }
    }
  }


  if (user->yavrg || user->zavrg)
  {
    printf("%d files used to calculate totals.", fileCnt);
    concOutput(user, tis, fileCnt);
  }

  else if (user->xflag || user->xpartflag)
  {
    // Make a list of active particles and sort.
    std::vector<int> sortedpId;
    sortedpId.reserve(user->partMap.size());
    for (auto & it : user->partMap)
      if (it.second.Flag() == 2)
        sortedpId.push_back(it.first);
    std::sort(sortedpId.begin(), sortedpId.end());

    if (user->xflag)
      crossXContOut(user, sortedpId);
    else if (user->xpartflag)
      crossXPartOut(user, sortedpId);
  }
  
  DAVecRestoreArray(user->fda, Count, &user->count);
  VecDestroy(Count);
  DAVecRestoreArray(user->fda, Coor, &user->coor);

  PetscFinalize();
}



PetscErrorCode ReadCoordinates(std::shared_ptr< UserCtx > & user)
{
  Cmpnts ***coor;

  Vec Coor;
  PetscInt bi, i, j, k, rank, IM, JM, KM;
  PetscReal *gc;
  FILE *fd;
  PetscReal       d0 = 1.;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscReal       cl = 1.;
  PetscOptionsGetReal(PETSC_NULL, "-cl", &cl, PETSC_NULL);

  char str[256];

  if(xyz_input) sprintf(str, "xyz.dat");
  else sprintf(str, "grid.dat");

  fd = fopen(str, "r");

  if(fd==NULL) printf("Cannot open %s !\n", str),exit(1);

  printf("Begin reading %s!\n", str);

  if(xyz_input) {i=1;}
  else if(binary_input) {
    fread(&i, sizeof(int), 1, fd);
    if(i!=1) PetscPrintf(PETSC_COMM_WORLD, "This seems to be a text file !\n"),exit(1);
  }
  else {
    fscanf(fd, "%i\n", &i);
    if(i!=1) PetscPrintf(PETSC_COMM_WORLD, "This seems to be a binary file !\n"),exit(1);
  }


  for (bi=block_number-1; bi>=0; bi--) {

    std::vector<double> X, Y,Z;
    double tmp;

    if(xyz_input) {
      fscanf(fd, "%i %i %i\n", &(/* user[bi]. */ user->IM), &(/* user[bi]. */ user->JM), &(/* user[bi]. */ user->KM));
      X.resize(/* user[bi]. */ user->IM);
      Y.resize(/* user[bi]. */ user->JM);
      Z.resize(/* user[bi]. */ user->KM);

      for (i=0; i</* user[bi]. */ user->IM; i++) fscanf(fd, "%le %le %le\n", &X[i], &tmp, &tmp);
      for (j=0; j</* user[bi]. */ user->JM; j++) fscanf(fd, "%le %le %le\n", &tmp, &Y[j], &tmp);
      for (k=0; k</* user[bi]. */ user->KM; k++) fscanf(fd, "%le %le %le\n", &tmp, &tmp, &Z[k]);
    }
    else if(binary_input) {
      fread(&(/* user[bi]. */ user->IM), sizeof(int), 1, fd);
      fread(&(/* user[bi]. */ user->JM), sizeof(int), 1, fd);
      fread(&(/* user[bi]. */ user->KM), sizeof(int), 1, fd);
    }
    else fscanf(fd, "%i %i %i\n", &(/* user[bi]. */ user->IM), &(/* user[bi]. */ user->JM), &(/* user[bi]. */ user->KM));

    IM = /* user[bi]. */ user->IM; JM = /* user[bi]. */ user->JM; KM = /* user[bi]. */ user->KM;


    DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_BOX,
               /* user[bi]. */ user->IM+1, /* user[bi]. */ user->JM+1, /* user[bi]. */ user->KM+1, 1,1,
               PETSC_DECIDE, 1, 2, PETSC_NULL, PETSC_NULL, PETSC_NULL,
               &(/* user[bi]. */ user->da));
    if(rans) {
      DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_BOX,
                 /* user[bi]. */ user->IM+1, /* user[bi]. */ user->JM+1, /* user[bi]. */ user->KM+1, 1,1,
                 PETSC_DECIDE, 2, 2, PETSC_NULL, PETSC_NULL, PETSC_NULL,
                 &(/* user[bi]. */ user->fda2));
    }
    DASetUniformCoordinates(/* user[bi]. */ user->da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    DAGetCoordinateDA(/* user[bi]. */ user->da, &(/* user[bi]. */ user->fda));

    DAGetLocalInfo(/* user[bi]. */ user->da, &(/* user[bi]. */ user->info));

    DALocalInfo     info = /* user[bi]. */ user->info;
    PetscInt        xs = info.xs, xe = info.xs + info.xm;
    PetscInt        ys = info.ys, ye = info.ys + info.ym;
    PetscInt        zs = info.zs, ze = info.zs + info.zm;
    PetscInt        mx = info.mx, my = info.my, mz = info.mz;

    DAGetGhostedCoordinates(/* user[bi]. */ user->da, &Coor);
    DAVecGetArray(/* user[bi]. */ user->fda, Coor, &coor);

    double buffer;

    for (k=0; k<KM; k++)
      for (j=0; j<JM; j++)
        for (i=0; i<IM; i++) {
          if(xyz_input) {}
          else if(binary_input) fread(&buffer, sizeof(double), 1, fd);
          else fscanf(fd, "%le", &buffer);

          if( k>=zs && k<=ze && j>=ys && j<ye && i>=xs && i<xe ) {
            if(xyz_input) coor[k][j][i].x = X[i]/cl;
            else coor[k][j][i].x = buffer/cl;
          }
        }

    for (k=0; k<KM; k++)
      for (j=0; j<JM; j++)
        for (i=0; i<IM; i++) {
          if(xyz_input) {}
          else if(binary_input) fread(&buffer, sizeof(double), 1, fd);
          else fscanf(fd, "%le", &buffer);

          if( k>=zs && k<=ze && j>=ys && j<ye && i>=xs && i<xe ) {
            if(xyz_input) coor[k][j][i].y = Y[j]/cl;
            else coor[k][j][i].y = buffer/cl;
          }
        }

    for (k=0; k<KM; k++)
      for (j=0; j<JM; j++)
        for (i=0; i<IM; i++) {
          if(xyz_input) {}
          else if(binary_input) fread(&buffer, sizeof(double), 1, fd);
          else fscanf(fd, "%le", &buffer);

          if( k>=zs && k<=ze && j>=ys && j<ye && i>=xs && i<xe ) {
            if(xyz_input) coor[k][j][i].z = Z[k]/cl;
            else coor[k][j][i].z = buffer/cl;
          }
        }

    DAVecRestoreArray(/* user[bi]. */ user->fda, Coor, &coor);

    Vec     gCoor;
    DAGetCoordinates(/* user[bi]. */ user->da, &gCoor);
    DALocalToGlobal(/* user[bi]. */ user->fda, Coor, INSERT_VALUES, gCoor);
    DAGlobalToLocalBegin(/* user[bi]. */ user->fda, gCoor, INSERT_VALUES, Coor);
    DAGlobalToLocalEnd(/* user[bi]. */ user->fda, gCoor, INSERT_VALUES, Coor);

  }

  fclose(fd);

  printf("Finish reading %s!\n", str);

  for (bi=0; bi<block_number; bi++) {
    /* user[bi]. */ user->_this = bi;
  }
  return(0);
}

void Calc_avg_shear_stress(std::shared_ptr< UserCtx > & user)
{
  double N=(double)tis+1.0;
  DALocalInfo   info = user->info;
  PetscInt      xs = info.xs, xe = info.xs + info.xm;
  PetscInt      ys = info.ys, ye = info.ys + info.ym;
  PetscInt      zs = info.zs, ze = info.zs + info.zm;
  PetscInt      mx = info.mx, my = info.my, mz = info.mz;

  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  PetscInt i, j, k;
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  Cmpnts ***usum, ***csi, ***eta, ***zet;
  PetscReal ***aj, ***psum, ***nvert;
  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;

  char filen[256];
  PetscViewer     viewer;

  Vec P_sum;
  DACreateGlobalVector(user->da, &P_sum);
  DACreateGlobalVector(user->fda, &user->Ucat_sum);

  ti=tis;
  sprintf(filen, "su0_%06d_%1d.dat", ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  VecLoadIntoVector(viewer, (user->Ucat_sum));
  PetscViewerDestroy(viewer);

  ti=tis;
  sprintf(filen, "sp_%06d_%1d.dat", ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  VecLoadIntoVector(viewer, P_sum);
  PetscViewerDestroy(viewer);

  DAVecGetArray(user->fda, user->Csi, &csi);
  DAVecGetArray(user->fda, user->Eta, &eta);
  DAVecGetArray(user->fda, user->Zet, &zet);
  DAVecGetArray(user->da, user->Aj, &aj);
  DAVecGetArray(user->da, user->Nvert, &nvert);
  DAVecGetArray(user->fda, user->Ucat_sum, &usum);
  DAVecGetArray(user->da, P_sum, &psum);


  double force_skin_bottom = 0;
  double force_pressure_bottom = 0;
  double force_bottom = 0;
  double area_bottom = 0;

  double force_skin_top = 0;
  double force_pressure_top = 0;
  double force_top = 0;
  double area_top = 0;

  j=0;
  for (k=lzs; k<lze; k++)
    for (i=lxs; i<lxe; i++) {
      if (nvert[k][j+1][i] < 0.1) {
        double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
        double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;

        dudc=0, dvdc=0, dwdc=0;

        dude=usum[k][j+1][i].x * 2.0 / N;
        dvde=usum[k][j+1][i].y * 2.0 / N;
        dwde=usum[k][j+1][i].z * 2.0 / N;

        dudz=0, dvdz=0, dwdz=0;

        double ajc = aj[k][j+1][i];
        double csi0 = csi[k][j+1][i].x, csi1 = csi[k][j+1][i].y, csi2 = csi[k][j+1][i].z;
        double eta0 = eta[k][j+1][i].x, eta1 = eta[k][j+1][i].y, eta2 = eta[k][j+1][i].z;
        double zet0 = zet[k][j+1][i].x, zet1 = zet[k][j+1][i].y, zet2 = zet[k][j+1][i].z;

        Compute_du_dxyz (csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc,
                         dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
                         &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );

        double j_area = sqrt( eta[k][j+1][i].x*eta[k][j+1][i].x + eta[k][j+1][i].y*eta[k][j+1][i].y + eta[k][j+1][i].z*eta[k][j+1][i].z );
        double ni[3], nj[3], nk[3];
        double nx, ny, nz;
        Calculate_normal(csi[k][j+1][i], eta[k][j+1][i], zet[k][j+1][i], ni, nj, nk);
        nx = nj[0]; //inward normal
        ny = nj[1]; //inward normal
        nz = nj[2]; //inward normal


        double Fp = - psum[k][j+1][i] * eta2 / N;
        double Fs = (dw_dx * nx + dw_dy * ny + dw_dz * nz) / user->ren * j_area;
        //double Fs = (du_dx * nx + du_dy * ny + du_dz * nz) / user->ren * j_area;

        force_skin_bottom += Fs;
        force_pressure_bottom += Fp;
        force_bottom += Fs + Fp;
        area_bottom += fabs(eta1);      // projected area
      }
    }

  j=my-2;
  for (k=lzs; k<lze; k++)
    for (i=lxs; i<lxe; i++) {
      if (nvert[k][j][i] < 0.1) {
        double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
        double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;

        dudc=0, dvdc=0, dwdc=0;

        dude = -usum[k][j][i].x * 2.0 / N;
        dvde = -usum[k][j][i].y * 2.0 / N;
        dwde = -usum[k][j][i].z * 2.0 / N;

        dudz=0, dvdz=0, dwdz=0;

        double ajc = aj[k][j][i];
        double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
        double eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
        double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;

        Compute_du_dxyz (csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc,
                         dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
                         &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );

        double j_area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
        double ni[3], nj[3], nk[3];
        double nx, ny, nz;
        Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
        nx = -nj[0]; //inward normal
        ny = -nj[1]; //inward normal
        nz = -nj[2]; //inward normal


        double Fp = - psum[k][j][i] * eta2 / N;
        double Fs = (dw_dx * nx + dw_dy * ny + dw_dz * nz) / user->ren * j_area;
        //double Fs = (du_dx * nx + du_dy * ny + du_dz * nz) / user->ren * j_area;

        force_skin_top += Fs;
        force_pressure_top += Fp;
        force_top += Fs + Fp;
        area_top += fabs(eta1); // projected area
      }
    }

  DAVecRestoreArray(user->fda, user->Csi, &csi);
  DAVecRestoreArray(user->fda, user->Eta, &eta);
  DAVecRestoreArray(user->fda, user->Zet, &zet);
  DAVecRestoreArray(user->da, user->Aj, &aj);
  DAVecRestoreArray(user->da, user->Nvert, &nvert);
  DAVecRestoreArray(user->fda, user->Ucat_sum, &usum);
  DAVecRestoreArray(user->da, P_sum, &psum);

  VecDestroy(P_sum);
  VecDestroy(user->Ucat_sum);

  printf("Top:\tarea=%f, force=%f, skin force=%f, pressure force=%f\n",
         area_top, force_top, force_skin_top, force_pressure_top);

  printf("\tstress=%f, skin stress=%f, pressure stress=%f, u*=%f, Re*=%f\n",
         force_top/area_top, force_skin_top/area_top, force_pressure_top/area_top,
         sqrt(fabs(force_top/area_top)), sqrt(fabs(force_top/area_top))*user->ren);

  printf("\n");

  printf("Bottom:\tarea=%f, force=%f, skin force=%f, pressure force=%f\n",
         area_bottom, force_bottom, force_skin_bottom, force_pressure_bottom);

  printf("\tstress=%f, skin stress=%f, pressure stress=%f, u*=%f, Re*=%f\n",
         force_bottom/area_bottom, force_skin_bottom/area_bottom, force_pressure_bottom/area_bottom,
         sqrt(fabs(force_bottom/area_bottom)), sqrt(fabs(force_bottom/area_bottom))*user->ren);
}

PetscErrorCode Lambda2(std::shared_ptr< UserCtx > & user)
{
  DALocalInfo   info = user->info;
  PetscInt      xs = info.xs, xe = info.xs + info.xm;
  PetscInt      ys = info.ys, ye = info.ys + info.ym;
  PetscInt      zs = info.zs, ze = info.zs + info.zm;
  PetscInt      mx = info.mx, my = info.my, mz = info.mz;

  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  PetscInt i, j, k;
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  Cmpnts ***ucat, ***csi, ***eta, ***zet;
  PetscReal ***aj, ***q, ***nvert;
  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;

  DAVecGetArray(user->fda, user->Ucat, &ucat);
  DAVecGetArray(user->fda, user->Csi, &csi);
  DAVecGetArray(user->fda, user->Eta, &eta);
  DAVecGetArray(user->fda, user->Zet, &zet);

  DAVecGetArray(user->da, user->Aj, &aj);
  DAVecGetArray(user->da, user->Nvert, &nvert);
  DAVecGetArray(user->da, user->P, &q);

  PetscReal uc, vc, wc, ue, ve, we, uz, vz, wz;

  PetscReal s11, s12, s13, s21, s22, s23, s31, s32, s33;
  PetscReal d11, d12, d13, d21, d22, d23, d31, d32, d33;

  PetscReal w11, w12, w13, w21, w22, w23, w31, w32, w33;
  //PetscReal so, wo;
  PetscReal csi1, csi2, csi3, eta1, eta2, eta3, zet1, zet2, zet3;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {

        double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
        double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
        double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
        double eta0= eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
        double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
        double ajc = aj[k][j][i];

        Compute_du_center (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
        Compute_du_dxyz (       csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
                                &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );

        double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
        double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),    Syz = 0.5*(dv_dz + dw_dy);
        double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);


        w11 = 0;
        w12 = 0.5*(du_dy - dv_dx);
        w13 = 0.5*(du_dz - dw_dx);
        w21 = -w12;
        w22 = 0.;
        w23 = 0.5*(dv_dz - dw_dy);
        w31 = -w13;
        w32 = -w23;
        w33 = 0.;


        double S[3][3], W[3][3], D[3][3];

        D[0][0] = du_dx, D[0][1] = du_dy, D[0][2] = du_dz;
        D[1][0] = dv_dx, D[1][1] = dv_dy, D[1][2] = dv_dz;
        D[2][0] = dw_dx, D[2][1] = dw_dy, D[2][2] = dw_dz;

        S[0][0] = Sxx;
        S[0][1] = Sxy;
        S[0][2] = Sxz;

        S[1][0] = Syx;
        S[1][1] = Syy;
        S[1][2] = Syz;

        S[2][0] = Szx;
        S[2][1] = Szy;
        S[2][2] = Szz;

        W[0][0] = w11;
        W[0][1] = w12;
        W[0][2] = w13;
        W[1][0] = w21;
        W[1][1] = w22;
        W[1][2] = w23;
        W[2][0] = w31;
        W[2][1] = w32;
        W[2][2] = w33;

        // lambda-2
        double A[3][3], V[3][3], d[3];

        for(int row=0; row<3; row++)
          for(int col=0; col<3; col++) A[row][col]=0;

        for(int row=0; row<3; row++)
          for(int col=0; col<3; col++) {
            A[row][col] += S[row][0] * S[0][col];
            A[row][col] += S[row][1] * S[1][col];
            A[row][col] += S[row][2] * S[2][col];
          }

        for(int row=0; row<3; row++)
          for(int col=0; col<3; col++) {
            A[row][col] += W[row][0] * W[0][col];
            A[row][col] += W[row][1] * W[1][col];
            A[row][col] += W[row][2] * W[2][col];
          }

        if(nvert[k][j][i]<0.1) {
          eigen_decomposition(A, V, d);
          q[k][j][i] = d[1];
        }
        else q[k][j][i] = 1000.0;
/*
// delta criterion
double DD[3][3];
for(int row=0; row<3; row++)
for(int col=0; col<3; col++) DD[row][col]=0;

for(int row=0; row<3; row++)
for(int col=0; col<3; col++) {
DD[row][col] += D[row][0] * D[0][col];
DD[row][col] += D[row][1] * D[1][col];
DD[row][col] += D[row][2] * D[2][col];
}
double tr_DD = DD[0][0] + DD[1][1] + DD[2][2];
double det_D = D[0][0]*(D[2][2]*D[1][1]-D[2][1]*D[1][2])-D[1][0]*(D[2][2]*D[0][1]-D[2][1]*D[0][2])+D[2][0]*(D[1][2]*D[0][1]-D[1][1]*D[0][2]);

//double Q = -0.5*tr_DD;

double SS=0, WW=0;
for(int row=0; row<3; row++)
for(int col=0; col<3; col++) {
SS+=S[row][col]*S[row][col];
WW+=W[row][col]*W[row][col];
}
double Q = 0.5*(WW - SS);

double R = - det_D;
if(nvert[k][j][i]<0.1) {
q[k][j][i] = pow( 0.5*R, 2. )  + pow( Q/3., 3.);
}
else q[k][j][i] = -10;
if(q[k][j][i]<0) q[k][j][i]=-10;
*/
      }
    }
  }

  DAVecRestoreArray(user->fda, user->Ucat, &ucat);
  DAVecRestoreArray(user->fda, user->Csi, &csi);
  DAVecRestoreArray(user->fda, user->Eta, &eta);
  DAVecRestoreArray(user->fda, user->Zet, &zet);

  DAVecRestoreArray(user->da, user->Aj, &aj);
  DAVecRestoreArray(user->da, user->Nvert, &nvert);
  DAVecRestoreArray(user->da, user->P, &q);

  return 0;
}

PetscErrorCode QCriteria(std::shared_ptr< UserCtx > & user)
{

  DALocalInfo   info = user->info;
  PetscInt      xs = info.xs, xe = info.xs + info.xm;
  PetscInt      ys = info.ys, ye = info.ys + info.ym;
  PetscInt      zs = info.zs, ze = info.zs + info.zm;
  PetscInt      mx = info.mx, my = info.my, mz = info.mz;

  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  PetscInt i, j, k;
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  Cmpnts ***ucat, ***csi, ***eta, ***zet;
  PetscReal ***aj, ***q, ***nvert;
  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;

  DAVecGetArray(user->fda, user->Ucat, &ucat);
  DAVecGetArray(user->fda, user->Csi, &csi);
  DAVecGetArray(user->fda, user->Eta, &eta);
  DAVecGetArray(user->fda, user->Zet, &zet);

  DAVecGetArray(user->da, user->Aj, &aj);
  DAVecGetArray(user->da, user->Nvert, &nvert);
  DAVecGetArray(user->da, user->P, &q);

  PetscReal uc, vc, wc, ue, ve, we, uz, vz, wz;

  PetscReal s11, s12, s13, s21, s22, s23, s31, s32, s33;
  PetscReal d11, d12, d13, d21, d22, d23, d31, d32, d33;

  PetscReal w11, w12, w13, w21, w22, w23, w31, w32, w33;
  PetscReal so, wo;
  PetscReal csi1, csi2, csi3, eta1, eta2, eta3, zet1, zet2, zet3;


  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
        double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
        double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
        double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
        double eta0= eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
        double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
        double ajc = aj[k][j][i];

        Compute_du_center (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
        Compute_du_dxyz (       csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
                                &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );

        double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
        double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),    Syz = 0.5*(dv_dz + dw_dy);
        double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
        so = Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz;

        w11 = 0;
        w12 = 0.5*(du_dy - dv_dx);
        w13 = 0.5*(du_dz - dw_dx);
        w21 = -w12;
        w22 = 0.;
        w23 = 0.5*(dv_dz - dw_dy);
        w31 = -w13;
        w32 = -w23;
        w33 = 0.;

        wo = w11*w11 + w12*w12 + w13*w13 + w21*w21 + w22*w22 + w23*w23 + w31*w31 + w32*w32 + w33*w33;

/*
  so = ( d11 *  d11 + d22 * d22 + d33 * d33) + 0.5* ( (d12 + d21) * (d12 + d21) + (d13 + d31) * (d13 + d31) + (d23 + d32) * (d23 + d32) );
  wo = 0.5 * ( (d12 - d21)*(d12 - d21) + (d13 - d31) * (d13 - d31) + (d23 - d32) * (d23 - d32) );
  V19=0.5 * ( (V13 - V11)*(V13 - V11) + (V16 - V12) * (V16 - V12) + (V17 - V15) * (V17 - V15) ) - 0.5 * ( V10 *  V10 + V14 * V14 + V18 * V18) - 0.25* ( (V13 + V11) * (V13 + V11) + (V16 + V12) * (V16 + V12) + (V17 + V15) * (V17 + V15) )
*/

        if( nvert[k][j][i]>0.1 ) q[k][j][i] = -100;
        else q[k][j][i] = (wo - so) / 2.;
      }
    }
  }

  DAVecRestoreArray(user->fda, user->Ucat, &ucat);
  DAVecRestoreArray(user->fda, user->Csi, &csi);
  DAVecRestoreArray(user->fda, user->Eta, &eta);
  DAVecRestoreArray(user->fda, user->Zet, &zet);

  DAVecRestoreArray(user->da, user->Aj, &aj);
  DAVecRestoreArray(user->da, user->Nvert, &nvert);
  DAVecRestoreArray(user->da, user->P, &q);

  return 0;
}

PetscErrorCode Velocity_Magnitude(std::shared_ptr< UserCtx > & user)        // store at P
{
  DALocalInfo   info = user->info;
  PetscInt      xs = info.xs, xe = info.xs + info.xm;
  PetscInt      ys = info.ys, ye = info.ys + info.ym;
  PetscInt      zs = info.zs, ze = info.zs + info.zm;
  PetscInt      mx = info.mx, my = info.my, mz = info.mz;

  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  PetscInt i, j, k;
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  Cmpnts ***ucat;
  PetscReal ***aj, ***q, ***nvert;
  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;

  DAVecGetArray(user->fda, user->Ucat, &ucat);
  DAVecGetArray(user->da, user->P, &q);

  PetscReal uc, vc, wc, ue, ve, we, uz, vz, wz;
  PetscReal d11, d12, d13, d21, d22, d23, d31, d32, d33;

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
        q[k][j][i] = sqrt( ucat[k][j][i].x*ucat[k][j][i].x + ucat[k][j][i].y*ucat[k][j][i].y + ucat[k][j][i].z*ucat[k][j][i].z );
      }
    }
  }

  DAVecRestoreArray(user->fda, user->Ucat, &ucat);
  DAVecRestoreArray(user->da, user->P, &q);

  return 0;
}



PetscErrorCode Velocity_Magnitude_4(std::shared_ptr< UserCtx > & user,
                                    PetscInt titi)       // store at P
{
  DALocalInfo   info = user->info;
  PetscInt      xs = info.xs, xe = info.xs + info.xm;
  PetscInt      ys = info.ys, ye = info.ys + info.ym;
  PetscInt      zs = info.zs, ze = info.zs + info.zm;
  PetscInt      mx = info.mx, my = info.my, mz = info.mz;

  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  PetscInt i, j, k;
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  Cmpnts ***ucat;
  PetscReal ***aj, ***q, ***nvert;
  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;

  DAVecGetArray(user->fda, user->Ucat, &ucat);
  DAVecGetArray(user->da, user->P, &q);

  PetscReal uc, vc, wc, ue, ve, we, uz, vz, wz;
  PetscReal d11, d12, d13, d21, d22, d23, d31, d32, d33;

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
        if (titi==1)q[k][j][i] = ucat[k][j][i].x;
        if (titi==2)q[k][j][i] = ucat[k][j][i].y;
        if (titi==3)q[k][j][i] = ucat[k][j][i].z;
      }
    }
  }

  DAVecRestoreArray(user->fda, user->Ucat, &ucat);
  DAVecRestoreArray(user->da, user->P, &q);

  return 0;
}






PetscErrorCode ibm_read(IBMNodes *ibm)
{
  PetscInt      rank;
  PetscInt      n_v , n_elmt ;
  PetscReal     *x_bp , *y_bp , *z_bp ;
  PetscInt      *nv1 , *nv2 , *nv3 ;
  PetscReal     *nf_x, *nf_y, *nf_z;
  PetscInt      i;
  PetscInt      n1e, n2e, n3e;
  PetscReal     dx12, dy12, dz12, dx13, dy13, dz13;
  PetscReal     t, dr;
  double xt;
  //MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(!rank) { // root processor read in the data
    FILE *fd;
    fd = fopen("ibmdata0", "r"); if (!fd) SETERRQ(1, "Cannot open IBM node file")
                                            n_v =0;
    fscanf(fd, "%i", &n_v);
    fscanf(fd, "%i", &n_v);
    fscanf(fd, "%le", &xt);
    ibm->n_v = n_v;

    MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    //    PetscPrintf(PETSC_COMM_WORLD, "nv, %d %e \n", n_v, xt);
    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));

    for (i=0; i<n_v; i++) {
      fscanf(fd, "%le %le %le %le %le %le", &x_bp[i], &y_bp[i], &z_bp[i], &t, &t, &t);
      x_bp[i] = x_bp[i] / 28.;
      y_bp[i] = y_bp[i] / 28.;
      z_bp[i] = z_bp[i] / 28.;
    }
    ibm->x_bp0 = x_bp; ibm->y_bp0 = y_bp; ibm->z_bp0 = z_bp;

    for (i=0; i<n_v; i++) {
      PetscReal temp;
      temp = ibm->y_bp0[i];
      ibm->y_bp0[i] = ibm->z_bp0[i];
      ibm->z_bp0[i] = -temp;
    }


    MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    fscanf(fd, "%i\n", &n_elmt);
    ibm->n_elmt = n_elmt;
    MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);
    for (i=0; i<n_elmt; i++) {

      fscanf(fd, "%i %i %i\n", nv1+i, nv2+i, nv3+i);
      nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;
      //      PetscPrintf(PETSC_COMM_WORLD, "I %d %d %d\n", nv1[i], nv2[i], nv3[i]);
    }
    ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

    fclose(fd);

    for (i=0; i<n_elmt; i++) {
      n1e = nv1[i]; n2e =nv2[i]; n3e = nv3[i];
      dx12 = x_bp[n2e] - x_bp[n1e];
      dy12 = y_bp[n2e] - y_bp[n1e];
      dz12 = z_bp[n2e] - z_bp[n1e];

      dx13 = x_bp[n3e] - x_bp[n1e];
      dy13 = y_bp[n3e] - y_bp[n1e];
      dz13 = z_bp[n3e] - z_bp[n1e];

      nf_x[i] = dy12 * dz13 - dz12 * dy13;
      nf_y[i] = -dx12 * dz13 + dz12 * dx13;
      nf_z[i] = dx12 * dy13 - dy12 * dx13;

      dr = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i] + nf_z[i]*nf_z[i]);

      nf_x[i] /=dr; nf_y[i]/=dr; nf_z[i]/=dr;


    }

    ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;

/*     for (i=0; i<n_elmt; i++) { */
/*       PetscPrintf(PETSC_COMM_WORLD, "%d %d %d %d\n", i, nv1[i], nv2[i], nv3[i]); */
/*     } */
    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }
  else if (rank) {
    MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_v = n_v;

    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
/*     ibm->x_bp0 = x_bp;  ibm->y_bp0 = y_bp; ibm->z_bp0 = z_bp; */

    MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_elmt = n_elmt;

    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

    ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;
    ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z;

    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

/*     MPI_Bcast(&(ibm->nv1), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nv2), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nv3), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */

/*     MPI_Bcast(&(ibm->nf_x), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nf_y), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nf_z), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
  }

/*   MPI_Barrier(PETSC_COMM_WORLD); */
  return(0);
}

PetscErrorCode ibm_read_ucd(IBMNodes *ibm)
{
  PetscInt      rank;
  PetscInt      n_v , n_elmt ;
  PetscReal     *x_bp , *y_bp , *z_bp ;
  PetscInt      *nv1 , *nv2 , *nv3 ;
  PetscReal     *nf_x, *nf_y, *nf_z;
  PetscInt      i;
  PetscInt      n1e, n2e, n3e;
  PetscReal     dx12, dy12, dz12, dx13, dy13, dz13;
  PetscReal     t, dr;
  PetscInt      temp;
  double xt;
  char string[128];
  //MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(!rank) { // root processor read in the data
    FILE *fd;
    fd = fopen("ibmdata", "r");
    if (!fd) SETERRQ(1, "Cannot open IBM node file")
               n_v =0;

    if (fd) {
      fgets(string, 128, fd);
      fgets(string, 128, fd);
      fgets(string, 128, fd);

      fscanf(fd, "%i %i %i %i %i\n", &n_v, &n_elmt, &temp, &temp, &temp);

      ibm->n_v = n_v;
      ibm->n_elmt = n_elmt;

      MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      //    PetscPrintf(PETSC_COMM_WORLD, "nv, %d %e \n", n_v, xt);
      PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
      PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
      PetscMalloc(n_v*sizeof(PetscReal), &z_bp);

      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));


      PetscReal cl = 1.0;//1./0.0254;

      PetscOptionsGetReal(PETSC_NULL, "-chact_leng_valve", &cl, PETSC_NULL);

      for (i=0; i<n_v; i++) {
        fscanf(fd, "%i %le %le %le", &temp, &x_bp[i], &y_bp[i], &z_bp[i]);
        x_bp[i] = x_bp[i] / cl;
        y_bp[i] = y_bp[i] / cl;
        z_bp[i] = z_bp[i] / cl;

        ibm->x_bp[i] = x_bp[i];
        ibm->y_bp[i] = y_bp[i];
        ibm->z_bp[i] = z_bp[i];

        ibm->u[i].x = 0.;
        ibm->u[i].y = 0.;
        ibm->u[i].z = 0.;
      }
      ibm->x_bp0 = x_bp; ibm->y_bp0 = y_bp; ibm->z_bp0 = z_bp;

      MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);



      MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

      PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);
      char str[20];
      for (i=0; i<n_elmt; i++) {

        fscanf(fd, "%i %i %s %i %i %i\n", &temp, &temp, str, nv1+i, nv2+i, nv3+i);
        nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;

      }
      ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

      fclose(fd);
    }
    for (i=0; i<n_elmt; i++) {
      n1e = nv1[i]; n2e =nv2[i]; n3e = nv3[i];
      dx12 = x_bp[n2e] - x_bp[n1e];
      dy12 = y_bp[n2e] - y_bp[n1e];
      dz12 = z_bp[n2e] - z_bp[n1e];

      dx13 = x_bp[n3e] - x_bp[n1e];
      dy13 = y_bp[n3e] - y_bp[n1e];
      dz13 = z_bp[n3e] - z_bp[n1e];

      nf_x[i] = dy12 * dz13 - dz12 * dy13;
      nf_y[i] = -dx12 * dz13 + dz12 * dx13;
      nf_z[i] = dx12 * dy13 - dy12 * dx13;

      dr = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i] + nf_z[i]*nf_z[i]);

      nf_x[i] /=dr; nf_y[i]/=dr; nf_z[i]/=dr;


    }

    ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;

    /*     for (i=0; i<n_elmt; i++) { */
    /*       PetscPrintf(PETSC_COMM_WORLD, "%d %d %d %d\n", i, nv1[i], nv2[i], nv3[i]); */
    /*     } */
    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }
  else if (rank) {
    MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_v = n_v;

    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

    PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));

/*     ibm->x_bp0 = x_bp;  ibm->y_bp0 = y_bp; ibm->z_bp0 = z_bp; */

    MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

    for (i=0; i<ibm->n_v; i++) {
      ibm->x_bp[i] = ibm->x_bp0[i];
      ibm->y_bp[i] = ibm->y_bp0[i];
      ibm->z_bp[i] = ibm->z_bp0[i];

      ibm->u[i].x = 0.;
      ibm->u[i].y = 0.;
      ibm->u[i].z = 0.;
    }
    MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_elmt = n_elmt;

    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

    ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;
    ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z;

    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

/*     MPI_Bcast(&(ibm->nv1), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nv2), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nv3), n_elmt, MPI_INTEGER, 0, PETSC_COMM_WORLD); */

/*     MPI_Bcast(&(ibm->nf_x), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nf_y), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
/*     MPI_Bcast(&(ibm->nf_z), n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD); */
  }

/*   MPI_Barrier(PETSC_COMM_WORLD); */
  return(0);
}

PetscErrorCode Combine_Elmt(IBMNodes *ibm, IBMNodes *ibm0, IBMNodes *ibm1)
{

  PetscInt i;

  ibm->n_v = ibm0->n_v + ibm1->n_v;
  ibm->n_elmt = ibm0->n_elmt + ibm1->n_elmt;

  PetscInt n_v = ibm->n_v, n_elmt = ibm->n_elmt;

  for (i=0; i<ibm0->n_v; i++) {
    ibm->x_bp[i] = ibm0->x_bp[i];
    ibm->y_bp[i] = ibm0->y_bp[i];
    ibm->z_bp[i] = ibm0->z_bp[i];

    ibm->u[i] = ibm0->u[i];
    ibm->uold[i] = ibm0->uold[i];
    //    ibm->u[i].x = 0.;
/*     PetscPrintf(PETSC_COMM_WORLD, "Vel %e %e %e\n", ibm->u[i].x, ibm->u[i].y, ibm->u[i].z); */
  }
  for (i=0; i<ibm0->n_elmt; i++) {
    ibm->nv1[i] = ibm0->nv1[i];
    ibm->nv2[i] = ibm0->nv2[i];
    ibm->nv3[i] = ibm0->nv3[i];

    ibm->nf_x[i] = ibm0->nf_x[i];
    ibm->nf_y[i] = ibm0->nf_y[i];
    ibm->nf_z[i] = ibm0->nf_z[i];
  }

  for (i=ibm0->n_v; i<n_v; i++) {
    ibm->x_bp[i] = ibm1->x_bp[i-ibm0->n_v];
    ibm->y_bp[i] = ibm1->y_bp[i-ibm0->n_v];
    ibm->z_bp[i] = ibm1->z_bp[i-ibm0->n_v];
    ibm->u[i].x = 0.;
    ibm->u[i].y = 0.;
    ibm->u[i].z = 0.;
  }

  for (i=ibm0->n_elmt; i<n_elmt; i++) {
    ibm->nv1[i] = ibm1->nv1[i-ibm0->n_elmt] + ibm0->n_v;
    ibm->nv2[i] = ibm1->nv2[i-ibm0->n_elmt] + ibm0->n_v;
    ibm->nv3[i] = ibm1->nv3[i-ibm0->n_elmt] + ibm0->n_v;

    ibm->nf_x[i] = ibm1->nf_x[i-ibm0->n_elmt];
    ibm->nf_y[i] = ibm1->nf_y[i-ibm0->n_elmt];
    ibm->nf_z[i] = ibm1->nf_z[i-ibm0->n_elmt];
  }

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    if (ti == ti) {
      FILE *f;
      char filen[80];
      sprintf(filen, "surface%06d.dat",ti);
      f = fopen(filen, "w");
      PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n");
      PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", n_v, 1670-96);
      for (i=0; i<n_v; i++) {

        /*    ibm->x_bp[i] = ibm->x_bp0[i];
              ibm->y_bp[i] = ibm->y_bp0[i];
              ibm->z_bp[i] = ibm->z_bp0[i] + z0;*/
        PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm->x_bp[i], ibm->y_bp[i], ibm->z_bp[i]);
      }
      for (i=96; i<n_elmt; i++) {
        if (fabs(ibm->nf_z[i]) > 0.5 ||
            (fabs(ibm->nf_z[i]) < 0.5 &&
             (ibm->x_bp[ibm->nv1[i]] * ibm->x_bp[ibm->nv1[i]] +
              ibm->y_bp[ibm->nv1[i]] * ibm->y_bp[ibm->nv1[i]]) < 0.44*0.44)) {
          PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
        }
      }
      fclose(f);

      sprintf(filen, "leaflet%06d.dat",ti);
      f = fopen(filen, "w");
      PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n");
      PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", n_v, 96);
      for (i=0; i<n_v; i++) {

        /*    ibm->x_bp[i] = ibm->x_bp0[i];
              ibm->y_bp[i] = ibm->y_bp0[i];
              ibm->z_bp[i] = ibm->z_bp0[i] + z0;*/
        PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm->x_bp[i], ibm->y_bp[i], ibm->z_bp[i]);
      }
      for (i=0; i<96; i++) {
        if (fabs(ibm->nf_z[i]) > 0.5 ||
            (fabs(ibm->nf_z[i]) < 0.5 &&
             (ibm->x_bp[ibm->nv1[i]] * ibm->x_bp[ibm->nv1[i]] +
              ibm->y_bp[ibm->nv1[i]] * ibm->y_bp[ibm->nv1[i]]) < 0.44*0.44)) {
          PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
        }
      }
      fclose(f);

    }
  }

  return 0;
}

PetscErrorCode Elmt_Move(IBMNodes *ibm, std::shared_ptr< UserCtx > & user)
{
  PetscInt n_v = ibm->n_v, n_elmt = ibm->n_elmt;

  PetscReal rcx = -0.122, rcz = -0.32, z0 = 4.52;
  rcx = -0.09450115; rcz = -0.3141615; z0 = 4.47;
  PetscReal dz;
  dz = -0.031;
  rcz = rcz-dz;
  rcx = rcx - dz * sin(10./180.*3.1415926);
  PetscReal temp;
  PetscInt i;

  PetscInt n1e, n2e, n3e;
  PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
  for (i=0; i<n_v; i++) {
    ibm->x_bp_o[i] = ibm->x_bp[i];
    ibm->y_bp_o[i] = ibm->y_bp[i];
    ibm->z_bp_o[i] = ibm->z_bp[i];
  }

  angle =-angle * 3.1415926/180.;
  //angle = 0;
  for (i=0; i<n_v/2; i++) {
    ibm->x_bp[i] = (ibm->x_bp0[i] -0.01- rcx) * cos(angle) - (ibm->z_bp0[i] - rcz) * sin(angle) + rcx;
    ibm->y_bp[i] = ibm->y_bp0[i];
    ibm->z_bp[i] = (ibm->x_bp0[i] -0.01- rcx) * sin(angle) + (ibm->z_bp0[i] - rcz) * cos(angle) + z0 + rcz;

  }
  rcx = -rcx;
  for (i=n_v/2; i<n_v; i++) {
    ibm->x_bp[i] = (ibm->x_bp0[i] +0.01- rcx) * cos(-angle) - (ibm->z_bp0[i] - rcz) * sin(-angle) + rcx;
    ibm->y_bp[i] = ibm->y_bp0[i];
    ibm->z_bp[i] = (ibm->x_bp0[i] +0.01- rcx) * sin(-angle) + (ibm->z_bp0[i] - rcz) * cos(-angle) + z0 + rcz;
  }

  /* Rotate 90 degree */
  for (i=0; i<n_v; i++) {
    temp = ibm->y_bp[i];
    ibm->y_bp[i] = ibm->x_bp[i];
    ibm->x_bp[i] = temp;
  }
  for (i=0; i<n_elmt; i++) {
    n1e = ibm->nv1[i]; n2e =ibm->nv2[i]; n3e =ibm->nv3[i];
    dx12 = ibm->x_bp[n2e] - ibm->x_bp[n1e];
    dy12 = ibm->y_bp[n2e] - ibm->y_bp[n1e];
    dz12 = ibm->z_bp[n2e] - ibm->z_bp[n1e];

    dx13 = ibm->x_bp[n3e] - ibm->x_bp[n1e];
    dy13 = ibm->y_bp[n3e] - ibm->y_bp[n1e];
    dz13 = ibm->z_bp[n3e] - ibm->z_bp[n1e];

    ibm->nf_x[i] = dy12 * dz13 - dz12 * dy13;
    ibm->nf_y[i] = -dx12 * dz13 + dz12 * dx13;
    ibm->nf_z[i] = dx12 * dy13 - dy12 * dx13;

    dr = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i] +
              ibm->nf_z[i]*ibm->nf_z[i]);

    ibm->nf_x[i] /=dr; ibm->nf_y[i]/=dr; ibm->nf_z[i]/=dr;

/*       PetscPrintf(PETSC_COMM_WORLD, "NFZ %d %d %d %d %e\n", i, ibm->nv1[i], ibm->nv2[i], ibm->nv3[i], ibm->nf_z[i]); */
    //      PetscPrintf(PETSC_COMM_WORLD, "%le %le %le %le %le %le\n", x_bp[n1e], y_bp[n1e], ibm->x_bp0[n1e], ibm->y_bp0[n1e], x_bp[n3e], y_bp[n3e]);

  }
  if (ti>0) {
    for (i=0; i<n_v; i++) {
      //      ibm->uold[i] = ibm->u[i];

      ibm->u[i].x = (ibm->x_bp[i] - ibm->x_bp_o[i]) / user->dt;
      ibm->u[i].y = (ibm->y_bp[i] - ibm->y_bp_o[i]) / user->dt;
      ibm->u[i].z = (ibm->z_bp[i] - ibm->z_bp_o[i]) / user->dt;
    }
  }
  else {
    for (i=0; i<n_v; i++) {
      ibm->u[i].x = 0.;
      ibm->u[i].y = 0.;
      ibm->u[i].z = 0.;
    }
  }
  return 0;
}

PetscErrorCode Elmt_Move1(IBMNodes *ibm, std::shared_ptr< UserCtx > & user)
{
  PetscInt n_v = ibm->n_v, n_elmt = ibm->n_elmt;

  PetscReal rcx = -0.122, rcz = -0.32, z0 = 4.52;
  rcx = -0.09450115; rcz = -0.3141615; z0 = 4.47;
  PetscReal dz;
  dz = -0.031;
  rcz = rcz-dz;
  rcx = rcx - dz * sin(10./180.*3.1415926);
  PetscReal temp;
  PetscInt i;

  PetscInt n1e, n2e, n3e;
  PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
  for (i=0; i<n_v; i++) {
    ibm->x_bp_o[i] = ibm->x_bp[i];
    ibm->y_bp_o[i] = ibm->y_bp[i];
    ibm->z_bp_o[i] = ibm->z_bp[i];
  }

  angle =-angle * 3.1415926/180.;
  //angle = 0;
  for (i=0; i<n_v/2; i++) {
    ibm->x_bp[i] = (ibm->x_bp0[i] -0.01- rcx) * cos(angle) - (ibm->z_bp0[i] - rcz) * sin(angle) + rcx;
    ibm->y_bp[i] = ibm->y_bp0[i];
    ibm->z_bp[i] = (ibm->x_bp0[i] -0.01- rcx) * sin(angle) + (ibm->z_bp0[i] - rcz) * cos(angle) + z0 + rcz;

  }
  rcx = -rcx;
  for (i=n_v/2; i<n_v; i++) {
    ibm->x_bp[i] = (ibm->x_bp0[i] +0.01- rcx) * cos(-angle) - (ibm->z_bp0[i] - rcz) * sin(-angle) + rcx;
    ibm->y_bp[i] = ibm->y_bp0[i];
    ibm->z_bp[i] = (ibm->x_bp0[i] +0.01- rcx) * sin(-angle) + (ibm->z_bp0[i] - rcz) * cos(-angle) + z0 + rcz;
  }

  /* Rotate 90 degree */
  for (i=0; i<n_v; i++) {
    temp = ibm->y_bp[i];
    ibm->y_bp[i] = ibm->x_bp[i];
    ibm->x_bp[i] = temp;
  }
  for (i=0; i<n_elmt; i++) {
    n1e = ibm->nv1[i]; n2e =ibm->nv2[i]; n3e =ibm->nv3[i];
    dx12 = ibm->x_bp[n2e] - ibm->x_bp[n1e];
    dy12 = ibm->y_bp[n2e] - ibm->y_bp[n1e];
    dz12 = ibm->z_bp[n2e] - ibm->z_bp[n1e];

    dx13 = ibm->x_bp[n3e] - ibm->x_bp[n1e];
    dy13 = ibm->y_bp[n3e] - ibm->y_bp[n1e];
    dz13 = ibm->z_bp[n3e] - ibm->z_bp[n1e];

    ibm->nf_x[i] = dy12 * dz13 - dz12 * dy13;
    ibm->nf_y[i] = -dx12 * dz13 + dz12 * dx13;
    ibm->nf_z[i] = dx12 * dy13 - dy12 * dx13;

    dr = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i] +
              ibm->nf_z[i]*ibm->nf_z[i]);

    ibm->nf_x[i] /=dr; ibm->nf_y[i]/=dr; ibm->nf_z[i]/=dr;

/*       PetscPrintf(PETSC_COMM_WORLD, "NFZ %d %d %d %d %e\n", i, ibm->nv1[i], ibm->nv2[i], ibm->nv3[i], ibm->nf_z[i]); */
    //      PetscPrintf(PETSC_COMM_WORLD, "%le %le %le %le %le %le\n", x_bp[n1e], y_bp[n1e], ibm->x_bp0[n1e], ibm->y_bp0[n1e], x_bp[n3e], y_bp[n3e]);

  }
  if (ti>0) {
    for (i=0; i<n_v; i++) {
      ibm->uold[i] = ibm->u[i];

      ibm->u[i].x = (ibm->x_bp[i] - ibm->x_bp_o[i]) / user->dt;
      ibm->u[i].y = (ibm->y_bp[i] - ibm->y_bp_o[i]) / user->dt;
      ibm->u[i].z = (ibm->z_bp[i] - ibm->z_bp_o[i]) / user->dt;
    }
  }
  else {
    for (i=0; i<n_v; i++) {
      ibm->u[i].x = 0.;
      ibm->u[i].y = 0.;
      ibm->u[i].z = 0.;
    }
  }
  return 0;
}


/*****************************************************************/
#ifdef MAX
#undef MAX
#endif

#define MAX(a, b) ((a)>(b)?(a):(b))

#define n 3

static double hypot2(double x, double y) {
  return sqrt(x*x+y*y);
}

// Symmetric Householder reduction to tridiagonal form.

static void tred2(double V[n][n], double d[n], double e[n]) {

//  This is derived from the Algol procedures tred2 by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

  for (int j = 0; j < n; j++) {
    d[j] = V[n-1][j];
  }

  // Householder reduction to tridiagonal form.

  for (int i = n-1; i > 0; i--) {

    // Scale to avoid under/overflow.

    double scale = 0.0;
    double h = 0.0;
    for (int k = 0; k < i; k++) {
      scale = scale + fabs(d[k]);
    }
    if (scale == 0.0) {
      e[i] = d[i-1];
      for (int j = 0; j < i; j++) {
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
        V[j][i] = 0.0;
      }
    } else {

      // Generate Householder vector.

      for (int k = 0; k < i; k++) {
        d[k] /= scale;
        h += d[k] * d[k];
      }
      double f = d[i-1];
      double g = sqrt(h);
      if (f > 0) {
        g = -g;
      }
      e[i] = scale * g;
      h = h - f * g;
      d[i-1] = f - g;
      for (int j = 0; j < i; j++) {
        e[j] = 0.0;
      }

      // Apply similarity transformation to remaining columns.

      for (int j = 0; j < i; j++) {
        f = d[j];
        V[j][i] = f;
        g = e[j] + V[j][j] * f;
        for (int k = j+1; k <= i-1; k++) {
          g += V[k][j] * d[k];
          e[k] += V[k][j] * f;
        }
        e[j] = g;
      }
      f = 0.0;
      for (int j = 0; j < i; j++) {
        e[j] /= h;
        f += e[j] * d[j];
      }
      double hh = f / (h + h);
      for (int j = 0; j < i; j++) {
        e[j] -= hh * d[j];
      }
      for (int j = 0; j < i; j++) {
        f = d[j];
        g = e[j];
        for (int k = j; k <= i-1; k++) {
          V[k][j] -= (f * e[k] + g * d[k]);
        }
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
      }
    }
    d[i] = h;
  }

  // Accumulate transformations.

  for (int i = 0; i < n-1; i++) {
    V[n-1][i] = V[i][i];
    V[i][i] = 1.0;
    double h = d[i+1];
    if (h != 0.0) {
      for (int k = 0; k <= i; k++) {
        d[k] = V[k][i+1] / h;
      }
      for (int j = 0; j <= i; j++) {
        double g = 0.0;
        for (int k = 0; k <= i; k++) {
          g += V[k][i+1] * V[k][j];
        }
        for (int k = 0; k <= i; k++) {
          V[k][j] -= g * d[k];
        }
      }
    }
    for (int k = 0; k <= i; k++) {
      V[k][i+1] = 0.0;
    }
  }
  for (int j = 0; j < n; j++) {
    d[j] = V[n-1][j];
    V[n-1][j] = 0.0;
  }
  V[n-1][n-1] = 1.0;
  e[0] = 0.0;
}

// Symmetric tridiagonal QL algorithm.

static void tql2(double V[n][n], double d[n], double e[n]) {

//  This is derived from the Algol procedures tql2, by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

  for (int i = 1; i < n; i++) {
    e[i-1] = e[i];
  }
  e[n-1] = 0.0;

  double f = 0.0;
  double tst1 = 0.0;
  double eps = pow(2.0,-52.0);
  for (int l = 0; l < n; l++) {

    // Find small subdiagonal element

    tst1 = MAX(tst1,fabs(d[l]) + fabs(e[l]));
    int m = l;
    while (m < n) {
      if (fabs(e[m]) <= eps*tst1) {
        break;
      }
      m++;
    }

    // If m == l, d[l] is an eigenvalue,
    // otherwise, iterate.

    if (m > l) {
      int iter = 0;
      do {
        iter = iter + 1;  // (Could check iteration count here.)

        // Compute implicit shift

        double g = d[l];
        double p = (d[l+1] - g) / (2.0 * e[l]);
        double r = hypot2(p,1.0);
        if (p < 0) {
          r = -r;
        }
        d[l] = e[l] / (p + r);
        d[l+1] = e[l] * (p + r);
        double dl1 = d[l+1];
        double h = g - d[l];
        for (int i = l+2; i < n; i++) {
          d[i] -= h;
        }
        f = f + h;

        // Implicit QL transformation.

        p = d[m];
        double c = 1.0;
        double c2 = c;
        double c3 = c;
        double el1 = e[l+1];
        double s = 0.0;
        double s2 = 0.0;
        for (int i = m-1; i >= l; i--) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = hypot2(p,e[i]);
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * (c * g + s * d[i]);

          // Accumulate transformation.

          for (int k = 0; k < n; k++) {
            h = V[k][i+1];
            V[k][i+1] = s * V[k][i] + c * h;
            V[k][i] = c * V[k][i] - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;

        // Check for convergence.

      } while (fabs(e[l]) > eps*tst1);
    }
    d[l] = d[l] + f;
    e[l] = 0.0;
  }

  // Sort eigenvalues and corresponding vectors.

  for (int i = 0; i < n-1; i++) {
    int k = i;
    double p = d[i];
    for (int j = i+1; j < n; j++) {
      if (d[j] < p) {
        k = j;
        p = d[j];
      }
    }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (int j = 0; j < n; j++) {
        p = V[j][i];
        V[j][i] = V[j][k];
        V[j][k] = p;
      }
    }
  }
}

void eigen_decomposition(double A[n][n], double V[n][n], double d[n]) {
  double e[n];
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      V[i][j] = A[i][j];
    }
  }
  tred2(V, d, e);
  tql2(V, d, e);
}


// Code added to handle particles.

void LOCAL::genRtree(std::shared_ptr< UserCtx > & user)
{
  PetscPrintf(PETSC_COMM_WORLD, "Starting grid genRtree\n");
  Util::Timer timer(Util::Timer::Wall);

  // Generate rtree for grid;
  new SpatialIndexGrid(user, user->spatialIndexGrid);
  user->rtreeGrid.reset(new Geom::Rtree(user->spatialIndexGrid));

  Geom::Subset set(Geom::Subset::Vector);
  set.resize(Id(user->IM, user->JM, user->KM));
  PetscInt xs, xe, ys, ye, zs, ze;
  user->spatialIndexGrid->indexRanges(xs, xe, ys, ye, zs, ze);
  for (int k=zs; k<ze-2; ++k)
  for (int j=ys; j<ye-2; ++j)
  for (int i=xs; i<xe-2; ++i)
  {
    set.add((Geom::Id) Id(i, j, k));
  }
  user->rtreeGrid->addSubset(set);

  std::stringstream mess;
  mess << timer;
  PetscPrintf(PETSC_COMM_WORLD, "%s\n", mess.str().c_str());
  PetscPrintf(PETSC_COMM_WORLD, "Ending grid genRtree\n");
}

std::shared_ptr<std::ifstream> openPartFile(int ti)
{
  char filename[128];
  sprintf(filename, "Particle%06d_0.dat", ti);
  std::shared_ptr<std::ifstream> fs = std::make_shared<std::ifstream>(filename);

  // If it didn't open, abort with an error.
  if (!fs->good())
  {
    printf("Couldn't open file %s\n", filename);
    exit(1);
  }

  // The first two line are headers. Read them before returning.
  std::string line;
  std::getline(*fs, line);
  std::getline(*fs, line);
    
  return fs;
}

int readPartFile(std::ifstream & fs, std::string & line,
                   std::array<int, 30> & start, std::array<int, 30> & end)
{
  std::getline(fs, line);

  if (fs.eof())
  {
    return 0;
  }

  if (fs.bad())
  {
    printf("Error reading the particle file.\n");
    exit(1);
  }

  int cnt = 0;
  start.fill(std::numeric_limits<int>::max());
  end.fill(std::numeric_limits<int>::max());

  std::size_t curr = 0;

  while (true)
  {
    // find the first non-white space character.
    for ( ; curr < line.size(); ++curr)
    {
      if (line[curr] == ' ' || line[curr] == '\t' || line[curr] == '\n')
        continue;
      break;
    }
    // if at the end of line, we're done here.
    if (curr >= line.size())
      break;

    start[cnt] = curr;

    // find the end of current field
    curr = line.find_first_of(" \t\n", curr);
    if (curr == std::string::npos)
    {
      end[cnt] = line.size();
      ++cnt;
      break;
    }
    end[cnt] = curr;
    ++cnt;
  }

  return 1;
}

int getIntField(int field, std::string line,
                std::array<int, 30> start, std::array<int, 30> end)
{
  int i;
  sscanf(line.substr(start[field], end[field]-start[field]).c_str(), "%d", &i);
  return i;
}

double getDoubleField(int field, std::string line,
                      std::array<int, 30> start, std::array<int, 30> end)
{
  double d;
  sscanf(line.substr(start[field], end[field]-start[field]).c_str(), "%lf", &d);
  return d;
}

void calcPartCount(std::shared_ptr< UserCtx > & user, int ti)
{
  std::string line;
  std::array<int, 30> start;
  std::array<int, 30> end;
  std::shared_ptr<std::ifstream> fs = openPartFile(ti);
  while (readPartFile(*fs, line, start, end))
  {
    int out = getIntField(10, line, start, end);

    if (user->end)
    {
      // If particle out (field indexed 10 from 0) status is 7 (out at boundary),
      // add 1 in proper entry.
      if (out != 7)
        continue;
    }
    else
    {
      // Want to count all active particles, out == 0.
      if (out != 0)
        continue;
    }

    double x = getDoubleField(0, line, start, end);
    double y = getDoubleField(1, line, start, end);
    double z = getDoubleField(2, line, start, end);

    double dist;
    double searchRad = 3.0 * user->spatialIndexGrid->getMaxDelta();
    auto result =
      user->rtreeGrid->getClosestEntity(Vector3d(x, y, z), &dist, searchRad);
    int i, j, k;
    Id(result).ijk(i, j, k);
    user->count[k][j][i] += 1;
  }
}

void calcXCrossing(std::shared_ptr< UserCtx > & user, int ti)
{
  std::string line;
  std::array<int, 30> start;
  std::array<int, 30> end;
  std::shared_ptr<std::ifstream> fs = openPartFile(ti);
  while (readPartFile(*fs, line, start, end))
  {
    // Get particle id and use it to lookup particle in partMap.
    int pId = getIntField(12, line, start, end);
    auto & part = user->partMap[pId];
  
    int out = getIntField(10, line, start, end);
    if (out != 0)
    {
      part.makeInactive();
      continue;
    }

    double x = getDoubleField(0, line, start, end);
    double y = getDoubleField(1, line, start, end);
    double z = getDoubleField(2, line, start, end);

    if (x >= GeomX::Tol(user->x))
      part.addMaxLoc(x, y, z);
  }
}

void CrossX::addMinLoc(double x, double y, double z)
{
  // if a min/max pair has been found or the particle is not active, do nothing.
  if (flag < 2)
  {
    flag = 1;
    xmin = x;
    ymin = y;
    zmin = z;
  }
}

void crossXPartOut(std::shared_ptr< UserCtx > & user,
                   std::vector<int> & sortedpId)
{
  char filename[128];
  sprintf(filename, "Particle.xsecpart.%+05.1f.dat", user->xpart);
  printf("Generating file %s\n", filename);

  auto f = fopen(filename, "w");
  if (!f)
  {
    printf("Couldn't open file %s.\n", filename);
    exit(1);
  }

  fprintf(f, "Variables =\tX,\tY,\tZ\tPID\n"
             "ZONE\tI = 1,\tJ = %d,\tDATAPACKING = POINT\n",
          sortedpId.size());

  for (auto & it : sortedpId)
  {
    double x;
    double y;
    double z;
    auto & part = user->partMap[it];
    part.Min(x, y, z);
    if (z == std::numeric_limits<double>::max())
    {
      part.Max(x, y, z);
      x = user->x;
    }
    else
      part.interpolate(user->x, x, y, z);
    fprintf(f, "%.5e\t%.5e\t%.5e\t%d\n", x, y, z, it);
  }

  fclose(f);
}


void crossXContOut(std::shared_ptr< UserCtx > & user,
                   std::vector<int> & sortedpId)
{
  PetscInt xs, xe, ys, ye, zs, ze;
  user->spatialIndexGrid->indexRanges(xs, xe, ys, ye, zs, ze);
  for (int k=zs; k<ze-2; k++)
  for (int j=ys; j<ye-2; j++)
  for (int i=xs; i<xe-2; i++)
    user->count[k][j][i] = 0.0;

  // Need k later on in the subroutine
  int k;
  bool first = true;
  for (auto & it : sortedpId)
  {
    double x;
    double y;
    double z;
    auto & part = user->partMap[it];
    part.Min(x, y, z);
    if (z == std::numeric_limits<double>::max())
    {
      part.Max(x, y, z);
      x = user->x;
    }
    else
      part.interpolate(user->x, x, y, z);

    double dist;
    double searchRad = 3.0 * user->spatialIndexGrid->getMaxDelta();
    auto result =
      user->rtreeGrid->getClosestEntity(Vector3d(x, y, z), &dist, searchRad);
    int i, j;
    Id(result).ijk(i, j, k);
    if (first)
    {
      first = false;
      printf("k value is %d\n", k);
    }
    user->count[k][j][i] += 1;
  }


  if (user->partfile == 1 || user->partfile == 2)
  {
    // Count the number of cells with a non-zero count.
    int totCnt = 0;

    for (int j=ys; j<ye-2; j++)
    for (int i=xs; i<xe-2; i++)
      if (user->count[k][j][i] != 0)
        ++totCnt;

    // Nothing to report
    if (totCnt == 0)
      return;

    printf("totCnt: %d\n", totCnt);
    
    char filename[128];
    sprintf(filename, "Particle.xseccont.%+05.1f.dat", user->x);
    printf("Generating file %s\n", filename);

    FILE *fd;
    fd = fopen(filename, "w");
    if (!fd)
    {
      printf("\nCould not open output file %s.\n\n", filename);
      exit(1);
    }

    /*
      Particle019980_0.dat
        Variables = X, Y, Z, CNT
        ZONE I = 1, J = 120000, DATAPACKING = POINT
        2.664000e+01 3.735525e+00 8.851376e+00 1.069666e+00
    */
    
    fprintf(fd, "Variables =\tX,\tY,\tZ,\tCNT\n");
    fprintf(fd, "ZONE\tI = 1,\tJ = %d,\tDATAPACKING = POINT\n", totCnt);

    for (int j=ys; j<ye-2; j++)
    for (int i=xs; i<xe-2; i++)
      if (user->count[k][j][i] != 0)
      {
        double y = 0.25 * (user->coor[k][j][i].y     + user->coor[k][j][i+1].y
                         + user->coor[k][j+1][i+1].y + user->coor[k][j+1][i].y);
        double z = 0.25 * (user->coor[k][j][i].z     + user->coor[k][j][i+1].z
                         + user->coor[k][j+1][i+1].z + user->coor[k][j+1][i].z);
        fprintf(fd, "%f\t%f\t%f\t%.0f\n", user->x, y, z, user->count[k][j][i]);
      }

    fclose(fd);
  }


  if (user->partfile == 0 || user->partfile == 2)
  {
    INTEGER4 I;
    char filename[128];
    sprintf(filename, "Result.xseccont.%+05.1f.plt", user->x);
    printf("\nGenerating file %s\n", filename);

    TECInitialize(user, "X Y Z CNT", filename);
    TECOutputGrid(user);
    TECOutputVector(user, user->count);
    I = TECEND100();
    if (I != 0)
    {
      printf("TECEND100 returned error %d\n", I);
      exit(1);
    }
  }
}

void CrossX::addMaxLoc(double x, double y, double z)
{
  // if a min/max pair has been found or the particle is not active, do nothing.
  if (flag < 2)
  {
    flag = 2;
    xmax = x;
    ymax = y;
    zmax = z;
  }
}

void CrossX::makeInactive()
{
  // If a min/max pair has been found, this entry can't be made inactive.
  if (flag != 2)
    flag = 3;
}

void CrossX::interpolate(double xTar, double & x, double & y, double & z)
{
  x = xTar;
  double m = (xTar - xmin) / (xmax - xmin);
  y = m*ymax + (1.0 - m)*ymin;
  z = m*zmax + (1.0 - m)*zmin;
}


void clearCount(std::shared_ptr< UserCtx > & user)
{
  PetscInt xs, xe, ys, ye, zs, ze;
  user->spatialIndexGrid->indexRanges(xs, xe, ys, ye, zs, ze);
  for (int k=zs; k<ze-2; k++)
  for (int j=ys; j<ye-2; j++)
  for (int i=xs; i<xe-2; i++)
    user->count[k][j][i] = 0.0;
}


void calcCount(std::shared_ptr< UserCtx > & user, int ti)
{
  static bool first = true;
  // Need to reset if we're not calculating file totals
  if (user->yconc || user->zconc)
    first = true;
  
  std::string line;
  std::array<int, 30> start;
  std::array<int, 30> end;
  std::shared_ptr<std::ifstream> fs = openPartFile(ti);

  char filename[128];
  sprintf(filename, "Particle%06d_0.dat", ti);
  printf("\nReading file %s\n", filename);

  while (readPartFile(*fs, line, start, end))
  {
    int out = getIntField(10, line, start, end);
    // Only include ground particles.
    if (user->zavrg && user->groundFlag)
    {
      if (out >= 3 and out <= 7)
        continue;
    }
    else if (out != 0)
      continue;

    double x = getDoubleField(0, line, start, end);
    double y = getDoubleField(1, line, start, end);
    if (user->yconc || user->yavrg)
      y = user->y;

    double z = getDoubleField(2, line, start, end);
    // Only include ground particles.
    if (user->zavrg && user->groundFlag)
      if (z > user->z)
        continue;

    if (user->zconc || user->zavrg)
      z = user->z;

    double dist;
    double searchRad = 3.0 * user->spatialIndexGrid->getMaxDelta();
    auto result =
      user->rtreeGrid->getClosestEntity(Vector3d(x, y, z), &dist, searchRad);
    int i, j, k;
    Id(result).ijk(i, j, k);
    if (first)
    {
      first = false;
      if (user->yconc || user->yavrg)
        printf("for y = %f, i = %d\n", y, i);
      if (user->zconc || user->zavrg)
        printf("for z = %f, j = %d\n", z, j);
    }
    user->count[k][j][i] += 1;
  }
}


void concOutput(std::shared_ptr< UserCtx > & user, int ti, int fileCnt)
{
  PetscInt xs, xe, ys, ye, zs, ze;
  user->spatialIndexGrid->indexRanges(xs, xe, ys, ye, zs, ze);

  char const * desc;
  if (user->yconc)
    desc = "concy";
  else if (user->yavrg)
    desc = "avrgy";
  else if (user->zconc)
    desc = "concz";
  else if (user->zavrg && user->ground)
    desc = "avrgz.ground";
  else if (user->zavrg)
    desc = "avrgz";

  if (user->partfile)
  {
    // Count the number of cells with a non-zero count.
    int totCnt = 0;

    for (int k=zs; k<ze-2; k++)
    for (int j=ys; j<ye-2; j++)
    for (int i=xs; i<xe-2; i++)
      if (user->count[k][j][i] != 0)
        ++totCnt;

    // Nothing to report
    if (totCnt == 0)
      return;

    printf("totCnt: %d\n", totCnt);
    
    char filename[128];
    sprintf(filename, "Particle%06d.%s.dat", ti, desc);
    printf("Generating file %s\n", filename);
    
    FILE *fd;
    fd = fopen(filename, "w");
    if (!fd)
    {
      printf("\nCould not open output file %s.\n\n", filename);
      exit(1);
    }

    /*
      Particle019980_0.dat
        Variables = X, Y, Z, CNT
        ZONE I = 1, J = 120000, DATAPACKING = POINT
        2.664000e+01 3.735525e+00 8.851376e+00 1.069666e+00
    */
    
    fprintf(fd, "Variables =\tX,\tY,\tZ,\tCNT\tAVG\tLOGAVG\n");
    fprintf(fd, "ZONE\tI = 1,\tJ = %d,\tDATAPACKING = POINT\n", totCnt);

    for (int k=zs; k<ze-2; k++)
    for (int j=ys; j<ye-2; j++)
    for (int i=xs; i<xe-2; i++)
      if (user->count[k][j][i] != 0)
      {
        double x, y, z;

        x = 0.25 * (user->coor[k][j][i].x     + user->coor[k][j][i+1].x
                  + user->coor[k][j+1][i+1].x + user->coor[k][j+1][i].x);

        // yconc or yavrg
        if (user->yconc || user->yavrg)
        {
          y = user->y;
          z = 0.25 * (user->coor[k][j][i].z     + user->coor[k][j][i+1].z
                    + user->coor[k][j+1][i+1].z + user->coor[k][j+1][i].z);
        }

        // zconc or zavrg
        else if (user->zconc || user->zavrg)
        {
          y = 0.25 * (user->coor[k][j][i].y     + user->coor[k][j][i+1].y
                    + user->coor[k][j+1][i+1].y + user->coor[k][j+1][i].y);
          z = user->z;
        }

        double avrg = user->count[k][j][i] / fileCnt;
        if (user->scaleFlag)
          avrg *= user->scale;
        double lavrg = log(avrg) / log(10.0);
        fprintf(fd, "%f\t%f\t%f\t%0f\t%12.6g\t%12.6g\n", x, y, z,
                user->count[k][j][i], avrg, lavrg);
      }

    fclose(fd);
  }
  else
  {
    INTEGER4 I;
    char filename[128];
    sprintf(filename, "Result%06d.%s.plt", ti, desc);
    printf("\nGenerating file %s\n", filename);

    // Scale the output.
    for (int k=zs; k<ze-2; k++)
    for (int j=ys; j<ye-2; j++)
    for (int i=xs; i<xe-2; i++)
      if (user->count[k][j][i] != 0)
      {
        user->count[k][j][i] /= fileCnt;
        if (user->scaleFlag)
          user->count[k][j][i] *= user->scale;
      }

    TECInitialize(user, "X Y Z CNT", filename);
    TECOutputGrid(user);
    TECOutputVector(user, user->count);
    I = TECEND100();
    if (I != 0)
    {
      printf("TECEND100 returned error %d\n", I);
      exit(1);
    }
  }
}


void calcZMaxCount(std::shared_ptr< UserCtx > & user, int ti)
{
  PetscInt xs, xe, ys, ye, zs, ze;
  user->spatialIndexGrid->indexRanges(xs, xe, ys, ye, zs, ze);
  for (int k=zs; k<ze-2; k++)
  for (int j=ys; j<ye-2; j++)
  for (int i=xs; i<xe-2; i++)
    user->count[k][j][i] = 0.0;

  bool first = true;

  std::string line;
  std::array<int, 30> start;
  std::array<int, 30> end;
  std::shared_ptr<std::ifstream> fs = openPartFile(ti);

  // j is needed below
  int j;
  
  while (readPartFile(*fs, line, start, end))
  {
    int out = getIntField(10, line, start, end);
    if (out != 0 && out != 8)
      continue;

    double x = getDoubleField(0, line, start, end);
    double y = getDoubleField(1, line, start, end);
    double z = getDoubleField(2, line, start, end);
    if (z > user->zmax)
      continue;
    z = user->zmax;

    double dist;
    double searchRad = 3.0 * user->spatialIndexGrid->getMaxDelta();
    auto result =
      user->rtreeGrid->getClosestEntity(Vector3d(x, y, z), &dist, searchRad);
    int i, k;
    Id(result).ijk(i, j, k);
    if (first)
    {
      first = false;
      printf("j value is %d\n", j);
    }
    user->count[k][j][i] += 1;
  }

  // Count the number of cells with a non-zero count.
  int totCnt = 0;

  for (int k=zs; k<ze-2; k++)
  for (int i=xs; i<xe-2; i++)
    if (user->count[k][j][i] != 0)
      ++totCnt;

  // Nothing to report
  if (totCnt == 0)
    return;

  printf("totCnt: %d\n", totCnt);
    
  char filename[128];
  sprintf(filename, "Particle%06d.zmaxcont.dat", ti);
  printf("Generating file %s\n", filename);

  FILE *fd;
  fd = fopen(filename, "w");
  if (!fd)
  {
    printf("\nCould not open output file %s.\n\n", filename);
    exit(1);
  }

  /*
    Particle019980_0.dat
      Variables = X, Y, Z, CNT
      ZONE I = 1, J = 120000, DATAPACKING = POINT
      2.664000e+01 3.735525e+00 8.851376e+00 1.069666e+00
  */
    
  fprintf(fd, "Variables =\tX,\tY,\tZ,\tCNT\n");
  fprintf(fd, "ZONE\tI = 1,\tJ = %d,\tDATAPACKING = POINT\n", totCnt);

  for (int k=zs; k<ze-2; k++)
  for (int i=xs; i<xe-2; i++)
    if (user->count[k][j][i] != 0)
    {
      double x = 0.25 * (user->coor[k][j][i].x     + user->coor[k][j][i+1].x
                       + user->coor[k][j+1][i+1].x + user->coor[k][j+1][i].x);
      double y = 0.25 * (user->coor[k][j][i].y     + user->coor[k][j][i+1].y
                       + user->coor[k][j+1][i+1].y + user->coor[k][j+1][i].y);
      fprintf(fd, "%f\t%f\t%f\t%.0f\n", x, y, user->z, user->count[k][j][i]);
    }

  fclose(fd);
}


void outputUandFTLE(std::shared_ptr< UserCtx > & user, int ti)
{
  INTEGER4 I;
  GridIndices idx(user);

  // check that ufield exists
  char ufilename[128];
  {
    sprintf(ufilename, "ufield%06d_0.dat", ti);
    auto f = fopen(ufilename, "r");
    if (f)
      fclose(f);
    else
    {
      printf("Couldn't open file %s.\n", ufilename);
      exit(1);
    }
  }
  
  // check that nvfield exists
  char const * nvfilename = "nvfield000000_0.dat";
  {

    // sprintf(nvfilename, "nvfield%06d_0.dat", ti);
    auto f = fopen(nvfilename, "r");
    if (f)
      fclose(f);
    else
    {
      printf("Couldn't open file %s.\n", nvfilename);
      exit(1);
    }
  }

  // Store ftle file names.
  std::vector < std::string > ftleFilenames;
  std::ifstream fstr("ftle.output");
  if (! fstr)
  {
    printf("\nCouldn't open file ftle.output so no ftle\n"
           "results will be included in this results file.\n");
  }
  else
  {
    std::string buff;
    while (std::getline(fstr, buff))
    {
      std::ifstream ftlef(buff);
      if (! ftlef)
        printf("\nCouldn't open file %s.\n"
               "It will not be included in the results file.", buff.c_str());
      else
        ftleFilenames.push_back(buff);
    }
    if (!ftleFilenames.size())
      printf("There where not valid filenames in the ftle.output file.\n"
             "No ftle data will be included in this results file/\n");
    fstr.close();
  }


  char filename[128];
  sprintf(filename, "Results%06d.uonly", ti);
  if (ftleFilenames.size())
    strcat(filename, ".ftle");
  strcat(filename, ".dat");
  printf("\nGenerating file %s\n", filename);

  
  auto header = std::string("X Y Z U V W UU NV");

  // Assumes the directory name has the follow format:
  // ftle[xyz].*9.9999
  for (auto const & it : ftleFilenames)
  {
    if (std::regex_search(it, std::regex("ftley")))
      header += " FTLEY";
    else if (std::regex_search(it, std::regex("ftlez")))
      header += " FTLEZ";
    else
      header += " FTLE";


    std::smatch cm;
    if (std::regex_search(it, cm, std::regex("([-0-9].[0-9]{4})")))
      header += cm[1];
    else
      header += &it-&ftleFilenames[0];
  }
  
  TECInitialize(user, header.c_str(), filename);

#if 0
  printf("Writing grid data.\n");
  TECOutputGrid(user);

  // output u,v,w
  // Everytime you create Ucat, it takes more memory so it's
  // important to recycle the memory.
  printf("Reading ufield data.\n");
  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, ufilename, FILE_MODE_READ, &viewer);
  if (!user->Ucat)
    DACreateGlobalVector(/*user[bi].*/ user->fda, &/*user[bi].*/ user->Ucat);
  VecLoadIntoVector(viewer, user->Ucat);
  PetscViewerDestroy(viewer);

  printf("Writing ufield data.\n");
  Cmpnts ***ucat;
  DAVecGetArray(/* user[bi]. */ user->fda, /* user[bi]. */ user->Ucat, &ucat);
  TECOutputCmpnts(idx.IMax, idx.JMax, idx.KMax,
                  idx.xs, idx.xe, idx.ys, idx.ye, idx.zs, idx.ze, ucat);
  DAVecRestoreArray(/*user[bi].*/ user->fda, /*user[bi].*/ user->Ucat, &ucat);

  // output nv
  // For most of my simulations, there is only one nvfield even
  // though they are generated for every timestep output.
  if (!user->Nvert)
  {
    printf("Reading nvfield data.\n");
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, nvfilename, FILE_MODE_READ, &viewer);
    DACreateGlobalVector(/*user[bi].*/ user->da, &/*user[bi].*/ user->Nvert);
    VecLoadIntoVector(viewer, user->Nvert);
    PetscViewerDestroy(viewer);
  }

  printf("Writing nvfield data.\n");
  PetscReal ***nvert;
  DAVecGetArray(/*user[bi].*/ user->da, /*user[bi].*/ user->Nvert, &nvert);
  TECOutputVector(user, nvert);
  DAVecRestoreArray(/*user[bi].*/ user->da, /*user[bi].*/ user->Nvert, &nvert);
#endif

  // output ftle files
  // Assume the following file format:
  // TITLE = FTLE FIELD
  // VARIABLES = X, Y, Z, FTLE
  // ZONE T='A150' I=157, J=833, F=POINT
  // -23.36000000  3.6865  -0.18000000  0.00000000
  for (auto const & it : ftleFilenames)
  {
    double dist;
    std::smatch cm;
    bool first = true;
    double searchRad = 3.0 * user->spatialIndexGrid->getMaxDelta();

    printf("Processing file %s.\n", it.c_str());

    char ftype = 'y';
    int axis[3];

    if (std::regex_search(it, std::regex("ftlex")))
    {
      ftype = 'x'; // k
      // axis = { 0, 1, 2 }
      axis[0] = 0;
      axis[1] = 1;
      axis[2] = 2;
    }
    else if (std::regex_search(it, std::regex("ftley")))
    {
      ftype = 'y'; // i
      // axis = { 1, 2, 0 }
      axis[0] = 1;
      axis[1] = 2;
      axis[2] = 0;
    }
    else if (std::regex_search(it, std::regex("ftlez")))
    {
      ftype = 'z'; // j
      // axis = { 2, 0, 1 }
      axis[0] = 2;
      axis[1] = 0;
      axis[2] = 1;
    }

    // Create the ftle array for tec output.
    // I'm going to store the ftle data in user->Count since it already exists.
    // It will need to be zeroed for the second or later timesteps so
    // we're going to just do it all the time.
    for (int k=idx.zs; k<idx.ze-2; k++)
    for (int j=idx.ys; j<idx.ye-2; j++)
    for (int i=idx.xs; i<idx.xe-2; i++)
      user->count[k][j][i] = 0.0;

    InterpolateVectorData ivd(idx.xs, idx.xe, idx.ys, idx.ye, idx.zs, idx.ze);
    
    std::string buff;
    std::ifstream fstr(it.c_str());
    while (std::getline(fstr, buff))
    {
      static std::regex redata(
        "^([-+]?[0-9]*.[0-9]+)[ \t]+([-+]?[0-9]*.[0-9]+)[ \t]+"
         "([-+]?[0-9]*.[0-9]+)[ \t]+([-+]?[0-9]*.[0-9]+)"
        );
      if (std::regex_search(buff, cm, redata))
      {
        double x = strtod(cm[1].str().c_str(), NULL);
        double y = strtod(cm[2].str().c_str(), NULL);
        double z = strtod(cm[3].str().c_str(), NULL);
        double ftle = strtod(cm[4].str().c_str(), NULL);
        
        auto result =
          user->rtreeGrid->getClosestEntity(Vector3d(x, y, z), &dist,searchRad);
        int ii, jj, kk;
        Id(result).ijk(ii, jj, kk);
        if (first)
        {
          first = false;
          if (ftype == 'x')
            printf("The k index corresponding to x=%f is %d.\n", x, kk);
          else if (ftype == 'y')
            printf("The i index corresponding to y=%f is %d.\n", y, ii);
          else if (ftype == 'z')
            printf("The j index corresponding to z=%f is %d.\n", z, jj);
        }

        // Update the values in count with this records ftle.
        ivd(axis, ii, jj, kk, Vector3d(x, y, z), user->coor, ftle, user->count);
      }
    }
    TECOutputVector(user, user->count);
  }
  
  I = TECEND100();
  if (I != 0)
  {
    printf("TECEND100 returned error %d\n", I);
    exit(1);
  }
}


// Interpolate data that is not exactly on the mesh onto the mesh.
// I'm using 1/d, where d is the distance of the known value to the
// mesh point, as the distribution weight.

// axis identifies which i,j,k indices are on the plane([0] and [1]) and
// which index is off the plane [2]. The calculations are done on the defined
// plane. Use 0 for i, 1 for j and 2 for k.

// i, j, k are the indices closeset to x,y,z
// x, y, z are the coordinates of the value location
// value is the value at x,y,z

InterpolateVectorData::InterpolateVectorData(int xs, int xe,
                                             int ys, int ye, int zs, int ze)
{
  // The start and end index values for the mesh in each direction.
  se[0][0] = xs;
  se[0][1] = xe-2;
  se[1][0] = ys;
  se[1][1] = ye-2;
  se[2][0] = zs;
  se[2][1] = ze-2;
}


void InterpolateVectorData::operator()(int const (&axis)[3],
                                       int const i, int const j, int const k,
                                       Vector3d const & c, Cmpnts ***coor,
                                       double val, PetscReal ***varr)
{
  double weightSum = 0.0;
  double weights[5][5] = { { 0.0, 0.0, 0.0, 0.0, 0.0 },
                           { 0.0, 0.0, 0.0, 0.0, 0.0 },
                           { 0.0, 0.0, 0.0, 0.0, 0.0 },
                           { 0.0, 0.0, 0.0, 0.0, 0.0 },
                           { 0.0, 0.0, 0.0, 0.0, 0.0 } };

  // Inverts transform os axis
  int iaxis[3];
  for (int ii=0; ii<3; ++ii)
    iaxis[axis[ii]] = ii;
  int const idx[3] = { i, j, k };

  // Calculate the value sums and weights.
  for (int ii=-2; ii<=2; ++ii)
  {
    printf("ii: %d\n", ii);
    // Check that is point is inside the mesh.
    if (idx[axis[0]]+ii <  se[axis[0]][0]) continue;
    if (idx[axis[0]]+ii >= se[axis[0]][1]) continue;
        
    for (int jj=-2; jj<=2; ++jj)
    {
      printf("jj: %d\n", jj);
      // Check that is point is inside the mesh.
      if (idx[axis[1]]+jj <  se[axis[1]][0]) continue;
      if (idx[axis[1]]+jj >= se[axis[1]][1]) continue;

      int it = idx[0] + (axis[0] == 0 ? ii : 0) + (axis[1] == 0 ? jj : 0);
      int jt = idx[1] + (axis[0] == 1 ? ii : 0) + (axis[1] == 1 ? jj : 0);
      int kt = idx[2] + (axis[0] == 2 ? ii : 0) + (axis[1] == 2 ? jj : 0);

      // When using axis values to get specific values, e.g. x, y or z, you need to use
      // the following indices:
      // for the x coordinate, use axis[2] (k)
      // for the y coordinate, use axis[0] (i)
      // for the z coordinate, use axis[1] (j)
      Vector3d coor2(coor[kt][jt][it]);
      double dist = sqrt(sqr(coor2(iaxis[0]) - c(iaxis[0]))
                       + sqr(coor2(iaxis[1]) - c(iaxis[1])));


      
      printf("axis:     %d,%d,%d\n", axis[0], axis[1], axis[2]);
      printf("ii,jj:    %d,%d\n", ii, jj);
      printf("i,j,k:    %d,%d,%d\n", idx[0], idx[1], idx[2]);
      printf("it,jt,kt: %d,%d,%d\n", it, jt, kt);
      printf("c:        %f,%f,%f\n", c(0), c(1), c(2));
      printf("coor2:    %f,%f,%f\n", coor2(0), coor2(1), coor2(2));
      

      
      // If the point is exactly on the mesh node, all the val goes there.
      if (dist == GeomX::Tol(0.0))
      {
        varr[kt][jt][it] += val;
        return;
      }

      weights[ii+2][jj+2] = 1.0 / dist;
      weightSum += weights[ii+2][jj+2];
    }
  }

  // Add the values to each mesh point.
  for (int ii=-2; ii<=2; ++ii)
  {
    // Check that is point is inside the mesh.
    if (idx[axis[0]]+ii <  se[axis[0]][0]) continue;
    if (idx[axis[0]]+ii >= se[axis[0]][1]) continue;

        
    for (int jj=-2; jj<=2; ++jj)
    {
      // Check that is point is inside the mesh.
      if (idx[axis[1]]+jj <  se[axis[1]][0]) continue;
      if (idx[axis[1]]+jj >= se[axis[1]][1]) continue;

      // Update this points value.
      if (weights[ii+2][jj+2] != 0.0)
      {
        int it = idx[0] + (axis[0] == 0 ? ii : 0) + (axis[1] == 0 ? jj : 0);
        int jt = idx[1] + (axis[0] == 1 ? ii : 0) + (axis[1] == 1 ? jj : 0);
        int kt = idx[2] + (axis[0] == 2 ? ii : 0) + (axis[1] == 2 ? jj : 0);
        varr[kt][jt][it] += val * weights[ii+2][jj+2] / weightSum;
      }
    }
  }
}
