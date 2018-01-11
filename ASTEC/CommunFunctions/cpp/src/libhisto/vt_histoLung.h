/*************************************************************************
 * vt_histoLung.h -
 *
 * $Id: vt_histoLung.h,v 1.2 2000/06/27 17:14:49 greg Exp $
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

#ifndef _vt_histoLung_h_
#define _vt_histoLung_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <vt_common.h>
#include <levenberg.h>


extern void _ProcessHistoLung( vt_image *theHisto,
			       double *theGaussParam );


#ifdef __cplusplus
}
#endif

#endif /* _vt_histoLung_h_  */
