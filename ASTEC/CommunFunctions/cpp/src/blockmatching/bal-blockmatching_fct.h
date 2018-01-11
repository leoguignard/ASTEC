/*************************************************************************
 * bal-blockmatching_fct.h -
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

#ifndef BAL_BLOCKMATCHING_FCT_H
#define BAL_BLOCKMATCHING_FCT_H

#include <bal-transformation-tools.h>
#include <bal-blockmatching-param.h>

#ifdef __cplusplus
extern "C" {
#endif

int blockmatching(
        char *floating_image,
        char *reference_image,
        char *result_image,
        char *initial_real_transformation,
        char *initial_voxel_transformation,
        char *initial_result_real_transformation,
        char *initial_result_voxel_transformation,
        char *result_real_transformation,
        char *result_voxel_transformation,
        int normalisation,
	bal_blockmatching_pyramidal_param param,
        int use_default_filename,
        char *command_line_file,
        char *log_file,
        int print_time,
        int print_parameters,
	int isDebug
        );

#ifdef __cplusplus
}
#endif

#endif //BAL_BLOCKMATCHING_FCT_H
