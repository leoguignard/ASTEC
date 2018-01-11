/*************************************************************************
 * vt_statsutil.h -
 *
 * $Id: vt_statsutil.h,v 1.1 2000/04/26 16:06:38 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Thu Mar 23 11:59:08 MET 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _vt_statsutil_h_
#define _vt_statsutil_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_image.h>
  
#include <vt_levenberg.h>
  
typedef struct {
  double sum0;
  double sum1;
  double sum2;
  double sum3;
  double sum4;
  double skewness;
  double kurtosis;
  double variance;
  double ecarttype;
  double moy;
  double med;
  int indMax;
  double valMax;

  double gauss[3];
  double nsgauss[4];
  
} typePoint;



#define _MINVAL_  1e-4


extern void _ComputeStatsValues( vt_image *imHisto,
				 typePoint *theXenon,
				 typePoint *theHmpao,
				 int minXenon, int maxXenon,
				 int minHmpao, int maxHmpao );


#ifdef __cplusplus
}
#endif

#endif /* _vt_statsutil_h_ */

