/*************************************************************************  
 * mt_membrane2D.h -
 *
 * $Id: mt_membrane2D.h,v 1.6 2013/06/20 11:28:34 gael Exp $
 *
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 * 
 * CREATION DATE: 
 * 2013/06/20
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _mt_membrane2D_h_
#define _mt_membrane2D_h_


#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
  /* #include <malloc.h> */

#include <vt_image.h>
#include <vt_common.h>
#include <recbuffer.h>

#include <vt_tubeutils.h>
#include <vt_tube2D.h>

/*
  vt_tube2D.h :
  typedef - mt_2Dtensor
          - vt_2Dimages
          - enumStructureColor
          - vt_3Dimres
          
  fonction- vt_2DtensorVoting
          - vt_2DtensorGaussianVoting
          - VT_Compute2DMultiScale
          - VT_Compute2DExtrema
          - VT_Compute2DMaskedExtrema
          - VT_Compute2DResponse
          - VT_FilterOn2DEigenValues    *
          - VT_Compute2DEigenVectors    *
          - VT_Filter2Dimages
          + VT_Alloc/Free/Write<...>    *
*/
typedef enum {
  TVCLASSIC,
  CFTV
} enumTVmode;

typedef enum {
  MEMBRANE,	//centerplanes
  VESSEL,	//centerlines
  BALL,		//centerballs
  MV,		//centerplanes+centerlines
  MB,		//centerplanes+centerballs
  VB,		//centerlines+centerballs
  MVB		//centerplanes+centerlines+centerballs
} shapeExtraction;

typedef enum {
  _REGULAR_,
  _RANDOM_
} MT_SAMPLINGMODE;


typedef enum {
  NONE,
  PLANE,
  LINE
} hessianMode;

typedef enum {
  SPARSE,
  DENSE
} enumVote ;

typedef struct {

  vt_image imxx;
  vt_image imyy;
  vt_image imxy;

  vt_image imvp1;
  vt_image imvp2;

  vt_image imtheta1;
  vt_image imtheta2;

  vt_image iszero;

} mt_2Dtensor;


extern int MT_Compute2DMultiScale( vt_image *theIm,
			    vt_3Dimres *imsRes,
			    double scale1,
			    double scale2,
			    int nbscales,
			    enumStructureColor color );

extern void MT_Compute2DExtrema( vt_3Dimres *imsRes,
			  vt_image *imExt );


extern void MT_Compute2DResponse( vt_2Dimages *par, enumStructureColor color,
			   double tau, double sigma );
			   
extern int MT_Filter2Dimages( vt_image *im, vt_2Dimages *par, float *theCoeffs);

extern int MT_SampleBin2D(vt_image *imageBin, double sample);

extern int MT_Compute2DTensorVoting(mt_2Dtensor *theTensor, vt_image **imagesIn, double scale, int niter,
        int nangles, enumTVmode TVmode, hessianMode initHessian,
        char *parName, int writeImages);

extern void MT_Compute2DTensorLineicExtrema( mt_2Dtensor *imTensor,
    vt_image *imExt );

extern void MT_Compute2DTensorBallExtrema( mt_2Dtensor *imTensor,
    vt_image *imExt );


extern void MT_Reinit2DTensorImage(mt_2Dtensor *theTensor);

extern void MT_Init2DTensorFromAngle(mt_2Dtensor *theTensor);

extern int MT_Compute2DStickFields( mt_2Dtensor **sfields, int *dimFields,
          double *angles, double *theCoeffs, int Nangles, enumTVmode mode);

extern int MT_Compute2DStickField( mt_2Dtensor *sfield, int *dimFields,
          double angle, double *theCoeffs, enumTVmode mode);

extern int MT_Compute2DBallFieldFromStick( mt_2Dtensor *bfield,
          mt_2Dtensor *sfields, int *dimFields,
          int Nangles);

extern int MT_Add2DField(mt_2Dtensor *tensorImg, mt_2Dtensor field, double coef,
          int x, int y, int slice, enumVote typeVote);

extern int set2DTensor(mt_2Dtensor *tensorImg, int x, int y, int slice,
          enumVote mode);

extern int nearest2DAngle(double theta, double *angles, int Nangles);

extern int MT_Compute2DAngles(double **angles, int Nangles);


extern int  VT_Alloc2DTensorFromImage( mt_2Dtensor *par, vt_image *im,
          char *genericname );
extern int  VT_Alloc2DTensor( mt_2Dtensor *par, char *n, int x, int y, int z,
          int t );
extern void VT_Free2DTensor ( mt_2Dtensor *par );
extern void VT_Write2DTensor( mt_2Dtensor *par );
extern void VT_Write2DTensorWithName( mt_2Dtensor *par, char *genericname);
extern int rand_(int a);

#ifdef __cplusplus
}
#endif

#endif /* _mt_membrane2D_h_ */
