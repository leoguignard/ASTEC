/*************************************************************************
 * vt_elfparam.h - extraction de parametres sur des parties numerotees
 *
 * $Id: vt_elfparam.h,v 1.10 2001/02/13 17:50:54 greg Exp $
 *
 * DESCRIPTION: 
 *
 *
 *
 *
 *
 * AUTHOR:
 * Gregoire Malandain
 *
 * CREATION DATE:
 * Jul 20 1999
 *
 * Copyright Gregoire Malandain, INRIA
 *
 *
 * ADDITIONS, CHANGES:
 *
 *
 */

#ifndef _vt_elfparam_h_
#define _vt_elfparam_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <vt_common.h>

extern void _VerboseInElfParam();
extern void _NoVerboseInElfParam();



typedef struct {

  int ptmin[3];
  int ptmax[3];

  int volume;
  int border;

  double equivSphereRadiusFromVolume;
  double equivSphereRadiusFromSurface;

  double maxDiameter;
  double medDiameter;
  double minDiameter;
  
  double maxDirection[3];
  double medDirection[3];
  double minDirection[3];

} typeComponentParameters;    






extern int _CreateArrayOfParametersFromImage( vt_image *theIm,
				       int slice,
				       typeComponentParameters **thePar );



extern void _ComputeComponentParameters( vt_image *theIm,
				  typeComponentParameters *thePar,
				  int slice,
				  int color,
				  int *corner,
				  int *windowSize );

extern void _GenerateEllipse2D( vt_image *theIm,
			 int slice,
			 int color,
			 double rmax,
			 double rmin,
			 double *radius,
			 double *angle );
extern void _GenerateEllipse3D( vt_image *theIm,
			 int color,
			 double rmax,
			 double rmin,
			 double *radius,
			 double *angle );



extern void _PrintEstimationErrors();
extern void _SetNbTestsForMaxFeretEstimation( int n );
extern void _SetFeretDiameterComputationToComputation();
extern void _SetFeretDiameterComputationToEstimation();
extern void _SetFeretDiameterComputationToComparison();



extern int _MaxValueInImage( vt_image *theIm,
		      int slice );
extern void _InitArrayOfParametersFromImage( vt_image *theIm,
				     typeComponentParameters *thePar,
				     int n );
extern void _FillArrayOfParametersFromImage( vt_image *theIm,
				      int slice,
				      typeComponentParameters *thePar );

extern void _TestsMaxDiameterInRandomList( int nbpoints,
					   int nbtests,
					   char *filename,
					   int testall );
#ifdef __cplusplus
}
#endif

#endif /* _vt_elfparam_h_ */
