/*************************************************************************
 * bal-transformation-list-tools.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mar 21 jan 2014 18:01:34 CET
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#ifndef BAL_TRANSFORMATION_LIST_TOOLS_H
#define BAL_TRANSFORMATION_LIST_TOOLS_H

#ifdef __cplusplus
extern "C" {
#endif





#include <bal-image.h>
#include <bal-transformation.h>


extern int BAL_ChangeTransformationList( bal_transformationList *theList,
					 bal_image *theIm,
					 bal_transformationList *resList,
					 bal_image *resIm,
					 int refTrsf,
					 int margin,
					 int iso );




#ifdef __cplusplus
}
#endif

#endif
