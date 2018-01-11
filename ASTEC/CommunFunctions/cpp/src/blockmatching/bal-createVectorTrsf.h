/*************************************************************************
 * applyTrsf.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2012, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mon Nov 19 17:45:00 CET 2012
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */

#ifndef BAL_CREATEVECTORTRSF_H
#define BAL_CREATEVECTORTRSF_H

#include <bal-transformation-tools.h>

#ifdef __cplusplus
extern "C" {
#endif

  typedef enum {
    SINUS_2D,
    SINUS_3D
  } enumVectorType;

  int createVectorTrsf(
		       char *restrsf_name,
		       char *template_image_name,
		       bal_integerPoint dim,
		       bal_doublePoint voxel,
		       enumTypeTransfo transformation_type,
		       enumVectorType vector_type,
		       int isDebug,
		       int isVerbose
		       );

#ifdef __cplusplus
}
#endif

#endif //BAL_CREATEVECTORTRSF_H
