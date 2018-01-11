/*************************************************************************
 * bal-composeTrsf.h -
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

#ifndef BAL_COMPOSETRSF_H
#define BAL_COMPOSETRSF_H

#include <bal-transformation-tools.h>

#ifdef __cplusplus
extern "C" {
#endif

  int composeTrsf(
		  char *restrsf_name,
		  char *template_image_name,
		  bal_integerPoint dim,
		  bal_doublePoint voxel,
		  /*
		    int first_trsf,
		    int last_trsf,
		  */
		  char** argv,
		  int argc,
		  int* is_a_trsf,
		  int _debug_,
		  int _verbose_
		  );

#ifdef __cplusplus
}
#endif

#endif //BAL_COMPOSETRSF_H
