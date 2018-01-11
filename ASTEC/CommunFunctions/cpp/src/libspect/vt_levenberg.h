/*************************************************************************
 * vt_levenberg.h -
 *
 * $Id: vt_levenberg.h,v 1.2 2000/04/20 08:07:04 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Wed Mar 22 09:23:04 MET 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _vt_levenberg_h_
#define _vt_levenberg_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <typedefs.h>
#include <systlin.h>


typedef double (*typeFuncEval)( double, double*, double* );	


extern double _LassenFunction( double x, double *thePar, double *derPar );
extern double _GaussianForLM ( double x, double *thePar, double *derPar );
extern double _NonSymetricGaussianForLM ( double x, double *thePar, double *derPar );



extern int VT_Modeling1DDataWithLevenberg( void *theX, bufferType xType,
				    void *theY, bufferType yType,
				    void *theW, bufferType wType,
				    void *theS, bufferType sType,
				    int theLength,
				    double *theParams, int nbParams,
				    typeFuncEval funcEval );


extern void VT_Levenberg_verbose ();
extern void VT_Levenberg_noverbose ();

#ifdef __cplusplus
}
#endif

#endif /* _vt_levenberg_h_  */
