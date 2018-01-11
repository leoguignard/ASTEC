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

#ifndef BAL_APPLYTRSF_H
#define BAL_APPLYTRSF_H

#include <bal-transformation-tools.h>

#ifdef __cplusplus
extern "C" {
#endif

  int applyTrsf(
		char *theim_name,
		char *resim_name,
		char *real_transformation_name,
		char *voxel_transformation_name,
		char *result_real_transformation_name,
		char *result_voxel_transformation_name,
		char *template_image_name,
		bal_integerPoint dim,
		bal_floatPoint voxel,
		int resize,
		enumTransformationInterpolation interpolation,
		ImageType type,
		int isDebug,
		int isVerbose
		);

#ifdef __cplusplus
}
#endif

#endif //BAL_APPLYTRSF_H
