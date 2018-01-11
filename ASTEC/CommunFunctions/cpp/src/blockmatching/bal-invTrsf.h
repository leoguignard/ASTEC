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

#ifndef BAL_INVTRSF_H
#define BAL_INVTRSF_H

#include <bal-transformation-tools.h>

#ifdef __cplusplus
extern "C" {
#endif

  int invTrsf(
	      char* thetrsf_name,
	      char* restrsf_name,
	      char* template_image_name,
	      bal_integerPoint dim,
	      bal_doublePoint voxel,
	      int isVerbose
	      );

#ifdef __cplusplus
}
#endif

#endif //BAL_INVTRSF_H
