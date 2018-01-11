/*****************************************************************************
 * vt_elfbarbules.h - enleve les barbules
 *
 * $Id: vt_elfbarbules.h,v 1.1 2000/02/29 10:43:38 greg Exp $
 *
 * Copyright©INRIA 2000
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * http://www.inria.fr/epidaure/personnel/malandain/
 * 
 * CREATION DATE: 
 * Feb, 24 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 *
 */

#ifndef _vt_elfbarbules_h_
#define _vt_elfbarbules_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>

#include <issimple3D.h>



extern int VT_RemoveBarbules( vt_image *theIm,
		       int labelMinCc,
		       int labelMaxCc,
		       int labelMinJc,
		       int labelMaxJc,
		       int length );




#ifdef __cplusplus
}
#endif

#endif 
