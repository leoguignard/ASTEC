/*************************************************************************
 * mt_membrane3D.h -
 *
 * $Id: mt_membrane3D.h,v 2.0 2013/10/22 11:10:00 gael Exp $
 *
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 * 
 * CREATION DATE: 
 * 2013/06/17
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _mt_membrane3D_h_
#define _mt_membrane3D_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <mt_membrane2D.h>
#include <vt_tube3D.h>

#define NELEMS(n) (sizeof(n) / sizeof (*n))


/*
  import vt_tube3D.h :
  typedef - typeResponseInSlice
          - typeCirclePoint
          - vt_3Dimages
          - vt_2Dimauxs
          - typeCSplinesInSlice
          
  fonction- VT_BuildCircle
          - VT_Compute3DresponseInSlice
          - VT_Filter3Dimages
          - VT_Compute3DMultiScale
          - VT_Compute3DExtrema
          - VT_Filter2Dimauxs
          - VT_ComputeNextRealFrame
          - VT_ComputeRealFrame
          - VT_ComputeRealFrameWithPivot
          - VT_ReechTriLinSlice
          - VT_ReechTriLinSlice2
          - VT_ReechCSplineSlice
          - VT_ReechCSplineSlice2
*/


typedef enum {
  MODELBASED,
  ACME
} enumMode;


typedef struct {

  vt_image imxx;
  vt_image imyy;
  vt_image imzz;
  vt_image imxy;
  vt_image imxz;
  vt_image imyz;

  vt_image imvp1;
  vt_image imvp2;
  vt_image imvp3;

  vt_image imtheta1;
  vt_image imphi1;
  vt_image imtheta2;
  vt_image imphi2;
  vt_image imtheta3;
  vt_image imphi3;
  
  vt_image iszero;

} vt_3Dtensor;


extern void MT_Compute3DresponseInSlice( typeResponseInSlice *aux,
					int dimx,
					int dimy,
					int dimz,
					int slice,
					double tau, float *theCoeff,
					 enumStructureColor color,
					 enumMode mode,
					 int hsp /* half size plane */);




extern int MT_Compute3DMultiScale( vt_image *theIm,
			    vt_3Dimres *imsRes,
			    double scale1,
			    double scale2,
			    int nbscales,
			    double zfact,
				   enumStructureColor color,
				   enumMode mode,
				   int hsp);

extern void MT_Compute3DExtrema( vt_3Dimres *imsRes,
                                 double zfact,
				 vt_image *imExt );


/* Calcul du gradient et des derivees selon z
    - calcul de x
    - calcul de y
    - calcul de z
    - calcul de zz
    - calcul de tmp0
    - calcul de tmp1
*/
extern int  MT_Filter3Dimages(vt_image *im, vt_3Dimages *par, float *theCoeffs);

/* Calcul de derivees dans le plan
   - calcul de xx
   - calcul de yy
   - calcul de xy
   - calcul de xz
   - calcul de yz
 */
extern int MT_Filter2Dimauxs( vt_3Dimages *ims, vt_2Dimauxs *par, 
			      float *theCoeffs, int slice );


extern int MT_SampleBin(vt_image *imageBin, double sample);

extern int MT_Compute3DTensorVoting( vt_3Dtensor *theTensor,
	  vt_image **imBin,  double scale, double zfact, int Niter, int Nangles, int Nsticks,
	  enumTVmode mode, hessianMode initHessian, char *parName, int writeImages );

extern int MT_Compute3DTensorVotingTest(vt_3Dtensor *theTensor, vt_image **imBin, int NanglesIter, double r,
        double alpha, hessianMode initHessian, char *parName, int writeImages);

extern int MT_Cumul3DTensorLine(vt_3Dtensor *theTensor, vt_image **imBin, double r, char *name);

extern int MT_AddLineToTensorImage3D(void *buf, 	// buffer array
                 int *dim, 		// buffer dimensions
                 bufferType t, 	// buffer type
                 int *pt, 		// origin point
                 double tht, 	// angle 1
                 double phi, 	// angle 2
                 int r,         // rayon de "vote"
                 double val     // value to be added
                );

extern int MT_ResetLineOfTensorImage3D(void *buf, int *dim, bufferType t, int *pt, double tht, double phi, int r);


extern void MT_ComputeTensorSurfaceExtrema( vt_3Dtensor *imTensor,
    vt_image *imExt, double zfact  );

extern void MT_ComputeTensorBallExtrema( vt_3Dtensor *imTensor,
    vt_image *imExt, double zfact );

extern void MT_ComputeTensorLineExtrema( vt_3Dtensor *imTensor,
    vt_image *imExt, double zfact );


extern void MT_ReinitTensorImage(vt_3Dtensor *theTensor);

extern void MT_InitTensorFromAngles(vt_3Dtensor *theTensor);

extern void MT_InitTensorTest(vt_3Dtensor *theTensor);

int MT_Compute3DTestFields( vt_3Dtensor **tfields, int *dimFields, 
          double **angles, int Nangles, double r, double alpha);

extern int MT_Compute3DTestField( vt_3Dtensor *tfield, int *dimFields, 
          double *angle, double ray, double Alpha);
                    
extern int MT_Compute3DStickFields( vt_3Dtensor **sfields, int *dimFields, 
          double **angles, double *theCoeffs, int Nangles, enumTVmode mode);

extern int MT_Compute3DStickField( vt_3Dtensor *sfield, int *dimFields, 
          double *angle, double *theCoeffs, enumTVmode mode);

extern int MT_Compute3DBallFieldFromStick( vt_3Dtensor *bfield,
          vt_3Dtensor *sfields, int *dimFields,
          int Nangles); 

extern int MT_Compute3DPlateFields( vt_3Dtensor **pfields, int *dimFields,
          double *theCoeffs, double **angles, int Nangles, int Nsticks,
          enumTVmode mode);

extern int MT_Compute3DPlateFieldFromSticks( vt_3Dtensor *pfield,
          vt_3Dtensor *sfields, int *dimFields, int Nsticks, int N);


extern int MT_AddField(vt_3Dtensor *tensorImg, vt_3Dtensor field, double coef,
          int x, int y, int z, enumVote typeVote);

extern int setTensor(vt_3Dtensor *tensorImg, int x, int y, int z,
          enumVote mode);

extern int nearestAngle(double theta, double phi, double **angles, int Nangles);


extern int  VT_Alloc3DtensorFromImage( vt_3Dtensor *par, vt_image *im,
          char *genericname );
extern int  VT_Alloc3Dtensor( vt_3Dtensor *par, char *n, int x, int y, int z,
          int t );
extern void VT_Free3Dtensor ( vt_3Dtensor *par );
extern void VT_Write3Dtensor( vt_3Dtensor *par );
extern void VT_Write3DtensorWithName( vt_3Dtensor *par, char *genericname);


typedef struct {
  double theta;
  double phi;

  double x;
  double y;
  double z;
} point;

typedef struct {
  point a;
  point b;
  point c;
} triangle;

extern void midPoint(point *mid, point a, point b);
extern void newTriangle(triangle *t, point a, point b, point c);

extern void setTriangle(triangle **T, triangle t, int n);
extern int addPoint(point **P, point p, int n, int n0);
extern int isequalPoint(point a, point b);
extern int MT_RecTriangulation(point **Pp, triangle **Tp, int *nT, int *nP,
          int niter);

extern int MT_Compute3DAngles(double ***angles, int Niter);

extern int MT_ComputeNormalAngles(double ***anglesPlane, int Nsticks,
          double theta, double phi);



extern int cmpAngles(double *angle1, double *angle2);



#ifdef __cplusplus
}
#endif

#endif /* _mt_membrane3D_h_ */

