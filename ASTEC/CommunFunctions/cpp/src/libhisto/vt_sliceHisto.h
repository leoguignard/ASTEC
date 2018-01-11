/*************************************************************************
 * vt_sliceHisto.h -
 *
 * $Id: vt_sliceHisto.h,v 1.9 2002/10/18 18:02:07 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Mon Jun 26 17:50:05 MET DST 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _vt_sliceHisto_h_
#define _vt_sliceHisto_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include <vt_common.h>




typedef enum {
  _AFFINE_,
  _LASSEN_
} enumIntensityCompensation;



typedef double (*typeIntensityTrsf)(double, double*);



int _GetReferenceSlice( int **theHisto, int dimz, int maxval );

double ** _GetPDFFromHisto( int **theHisto, int n, int max );

int ** _GetSlicesHisto( vt_image *theIm, int *max );
int * _GetImageHisto( vt_image *theIm, vt_image *immask, int *max );

void _PrintOneHistoForMatlab( int fd, FILE *fp, 
			      double *theCoeff, int *theHisto,
			      int length, int id, 
			      enumIntensityCompensation c ) ;


double _affine( double x, double *c );
double _inv_affine( double x, double *c );

double _constant( double x, double *c );
double _inv_constant( double x, double *c );


double _DirecteSSD( double *c, void *par );
double _InverseSSD( double *c, void *par );
double _SymetrieSSD( double *c, void *par );

double _DirecteCorrelation( double *c, void *par );
double _InverseCorrelation( double *c, void *par );
double _SymetrieCorrelation( double *c, void *par );

double _DirecteVraisemblance( double *c, void *par );
double _InverseVraisemblance( double *c, void *par );
double _SymetrieVraisemblance( double *c, void *par );
double _MinimumVraisemblance( double *c, void *par );

void _initAffineTrsfBetweenTwoHisto( double *theCoeff,
			       int *theHisto, 
			       int hlength,
			       int *refHisto,
			       int rlength );

void _initConstantTrsfBetweenTwoHisto( double *theCoeff,
			       int *theHisto, 
			       int hlength,
			       int *refHisto,
			       int rlength );

void _evalTrsfBetweenTwoHisto( double *theCoeff,
			       int *theHisto, 
			       int hlength,
			       int *refHisto,
			       int rlength,
			       typeIntensityTrsf trsf,
			       typeIntensityTrsf invtrsf,
			       int n,
			       double (*func)(double *, void *),
			       double sigma );

void _evalTrsfsIn3DHistos( double **theCoeff,
			   int **theHisto, 
			   int length,
			   int dimz,
			   int iref,
			   typeIntensityTrsf trsf,
			   typeIntensityTrsf invtrsf,
			   double (*func)(double *, void *),
			   double sigma );

void _multiScaleEvalTrsfBetweenTwoHisto( double *theCoeff,
					 int *theHisto, 
					 int hlength,
					 int *refHisto,
					 int rlength,
					 typeIntensityTrsf trsf,
					 typeIntensityTrsf invtrsf,
					 int nparam,
					 double (*func)(double *, void *),
					 double sigma,
					 int lscale, int nscale );








extern void _FitGaussiansOnSliceHisto( double *par,
				   int **theHisto, int dimz, int maxval,
				   int imin, int imax,
				   int z, int ng );

extern void _ComputeBiasFromJointHisto( vt_image *image, double ** theCoeff,
				 int iref );


extern void _ComputeBiasFromGaussians( FILE *fp,
				       double ** theCoeff,
				       int **theHisto, int dimz, int maxval,
				       int imin, int imax, int iref, int id );



extern double ** _GetPDFFromHisto( int **theHisto, int n, int max );
extern int ** _GetSlicesHisto( vt_image *theIm, int *max );



extern double _CorrelationOfOneHisto( double *c, void *par );
extern double _SSDOfOneHisto( double *c, void *par );
extern double _SADOfOneHisto( double *c, void *par );
extern double _entropyOfOneHisto( double *c, void *par );
extern void _minimizeFunctionWRTReference( double **theCoeff, int iref,
				    int **theHisto, int dimz, int maxval,
				    double (*func)(double *, void *) );



extern void _Print2DSlicesHistoForMatlab( int fd, FILE *fp, int **theHisto, int dimz, int maxval, int id );
extern void _Print2DSlicesPDFForMatlab( int fd, FILE *fp, double **theHisto, int dimz, int maxval, int id );
extern void _PrintOneSliceHistoForMatlab( int fd, FILE *fp, 
					  double ** theCoeff, int **theHisto, 
					  int slice, int maxval, int id );
extern void _PrintOneSlicePDFForMatlab( int fd, FILE *fp, 
					double ** theCoeff, double **theHisto, 
					int slice, int maxval, int id );


extern void _PrintJointHistoOfTwoSlicesForMatlab( int fd, FILE *fp, 
						  vt_image *image, int maxval,
						  int s1, int s2 );

extern void _PrintSlicesHistoForMatlab( int fd, FILE *fp, double ** theCoeff,
					int **theHisto, int dimz, int maxval, int id );




extern void _PrintFunctionOfTwoSlicesForMatlab( int fd, FILE *fp,
					int **theHisto, int dimz, int maxval,
					 double (*func)(double *, void *),
					int s1, int s2 );



typedef struct {
  double x;
  double xmin;
  double xmax;
  double p;
  
} typeProbabilite;


typeProbabilite * _buildOneNewPDF( double binsize, int maxval, int *nb );
void _newPDFFromPDF( typeProbabilite *pdf, int npdf,
		    double *thePdf, int n, double sigma );
void _PrintOnePDFForMatlab( int fd, FILE *fp, typeProbabilite *pdf, int n, int id );



#ifdef __cplusplus
}
#endif

#endif /* _vt_sliceHisto_h_  */
