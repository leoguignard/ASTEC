/*************************************************************************
 * vt_jointhisto.h -
 *
 * $Id: vt_jointhisto.h,v 1.2 2000/03/23 15:19:14 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Wed Mar 22 22:25:15 MET 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _vt_jointhisto_h_
#define _vt_jointhisto_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <math.h>

#include <vt_image.h>

extern int ComputeJointHistoWithTrsfAndMask( vt_image *imageHmpao,
				      vt_image *imageXenon,
				      vt_image *maskHmpao,
				      vt_image *maskXenon,
				      vt_image *histo,
				      double *mat,
				      float sigmaHmpao,
				      float sigmaXenon,
				      int minHmpao,
				      int minXenon );

extern int ComputeJointHistoWithTrsf( vt_image *imageHmpao,
				      vt_image *imageXenon,
				      vt_image *histo,
				      double *mat,
				      float sigmaHmpao,
				      float sigmaXenon,
				      int minHmpao,
				      int minXenon );

extern int ComputeJointHistoWithoutTrsf( vt_image *imageHmpao,
					 vt_image *imageXenon,
					 vt_image *histo,
					 float sigmaHmpao,
					 float sigmaXenon,
					 int minHmpao );


#ifdef __cplusplus
}
#endif

#endif /* _vt_jointhisto_h_ */

