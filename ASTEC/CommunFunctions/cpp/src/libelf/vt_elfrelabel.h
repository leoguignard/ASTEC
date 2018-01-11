/*****************************************************************************
 * vt_elfrelabel.h - renumerote les composantes connexes par taille croissante
 *
 * $Id: vt_elfrelabel.h,v 1.1 2000/02/29 10:43:39 greg Exp $
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

#ifndef _vt_elfrelabel_h_
#define _vt_elfrelabel_h_

#ifdef __cplusplus
extern "C" {
#endif


#include <typedefs.h>



extern int RelabelConnectedComponentsSortBySize( void *inputBuf,
						 bufferType typeIn,
						 int *theDim );




#ifdef __cplusplus
}
#endif

#endif 
