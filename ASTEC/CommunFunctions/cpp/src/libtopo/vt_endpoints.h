/*************************************************************************
 * vt_endpoints.h -
 *
 * $Id: vt_endpoints.h,v 1.2 2001/04/03 10:51:59 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Tue Mar 20 16:41:37 MET 2001
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _vt_endpoints_h_
#define _vt_endpoints_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <string.h>
#include <vt_common.h>




extern int VT_ComputeBackDistanceInsideObjects( vt_image *theCC,
						vt_image *theDist,
						vt_image *theBackDist,
						int *inc );


extern int VT_InitialiseImageToBeThinned( vt_image *theBackDist,
					  vt_image *theImToBeThinned,
					  int halfWindowSize,
					  int *inc );


#ifdef __cplusplus
}
#endif

#endif /* _vt_endpoints_h_ */

