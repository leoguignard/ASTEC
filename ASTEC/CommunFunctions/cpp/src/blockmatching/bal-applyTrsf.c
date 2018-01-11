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
#include "bal-applyTrsf.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <bal-stddef.h>
#include <bal-image.h>
#include <bal-transformation.h>
#include <bal-transformation-tools.h>

static char *program = "ApplyTrsf";

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
	      )
{
  /* Initialize parameters
   */
  bal_image theim;
  bal_image resim;
  bal_image template;
  bal_transformation theTrsf;

  bal_image auxim;

  BAL_InitImage( &theim, NULL, 0, 0, 0, 0, UCHAR );
  BAL_InitImage( &resim, NULL, 0, 0, 0, 0, UCHAR );
  BAL_InitImage( &template, NULL, 0, 0, 0, 0, UCHAR );
  BAL_InitTransformation( &theTrsf );


  /* reading input image
   */
  if ( BAL_ReadImage( &theim, theim_name, 0 ) != 1 ) {
    if ( isVerbose )
      fprintf( stderr, "%s: can not read input image '%s'\n", program, theim_name );
    return -1;
  }

  /* reading transformation, if any 
   */
  if ( real_transformation_name != NULL ) {
    if ( BAL_ReadTransformation( &theTrsf, real_transformation_name ) != 1 ) {
      BAL_FreeImage( &theim );
      if ( isVerbose )
	fprintf( stderr, "%s: unable to read 'real' transformation '%s'\n", program, real_transformation_name );
      return -1;
    }
  }
  else if ( voxel_transformation_name != NULL ) {
    if ( BAL_ReadTransformation( &theTrsf, voxel_transformation_name ) != 1 ) {
      BAL_FreeImage( &theim );
      if ( isVerbose )
	fprintf( stderr, "%s: unable to read 'voxel' transformation '%s'\n", program, voxel_transformation_name );
      return -1;
    }
    theTrsf.transformation_unit = VOXEL_UNIT;
  }



  /* initializing result image
     - with reference image, if any
     - with parameters, if any
     - with transformation, if vector field
     - with input image
  */

  /* initialisation with a reference image 
   */
  if ( template_image_name != NULL ) {
    if ( BAL_ReadImage( &template, template_image_name, 0 ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      BAL_FreeImage( &theim );
      if ( isVerbose )
	fprintf( stderr, "%s: can not read template image '%s'\n", program, template_image_name );
      return -1;
    }
    if ( BAL_InitAllocImage( &resim, (char*)NULL, 
			     template.ncols, template.nrows, template.nplanes, 
			     theim.vdim, theim.type ) != 1 ) {
      BAL_FreeImage( &template );
      BAL_FreeTransformation( &theTrsf );
      BAL_FreeImage( &theim );
      if ( isVerbose )
	fprintf( stderr, "%s: unable to initialize result image (from template image)\n", program );
      return -1;
    }
    resim.vx = template.vx;
    resim.vy = template.vy;
    resim.vz = template.vz;
    BAL_FreeImage( &template );
  }

  /* initialisation with parameters (dimensions)
   */
  else if ( dim.x > 0 && dim.y > 0 ) {
    if ( dim.z > 0 ) {
      if ( BAL_InitAllocImage( &resim, (char*)NULL, dim.x, dim.y, dim.z, 1, theim.type ) != 1 ) {
	BAL_FreeTransformation( &theTrsf );
	BAL_FreeImage( &theim );
	if ( isVerbose )
	  fprintf( stderr, "%s: unable to initialize result image (from parameters, 3D case)\n", program );
	return -1;
      }
    }
    else {
      if ( BAL_InitAllocImage( &resim, (char*)NULL, dim.x, dim.y, 1, 1, theim.type ) != 1 ) {
	BAL_FreeTransformation( &theTrsf );
	BAL_FreeImage( &theim );
	if ( isVerbose )
	  fprintf( stderr, "%s: unable to initialize result image (from parameters, 2D case)\n", program );
	return -1;
      }
    }
    if ( voxel.x > 0.0 ) resim.vx = voxel.x;
    if ( voxel.y > 0.0 ) resim.vy = voxel.y;
    if ( voxel.z > 0.0 ) resim.vz = voxel.z;
  }

  /* initialisation with parameters (voxel size)
   */
  else if ( resize && (voxel.x > 0 && voxel.y > 0) ) {
    dim.x = theim.ncols * theim.vx / voxel.x;
    dim.y = theim.nrows * theim.vy / voxel.y;
    
    if ( theim.nplanes > 1 ) {
      if ( voxel.z <= 0.0 ) {
	BAL_FreeTransformation( &theTrsf );
	BAL_FreeImage( &theim );
	if ( isVerbose )
	  fprintf( stderr, "%s: unable to calculate result image Z dimension (from parameters, 3D case)\n", program );
	return -1;
      }
      dim.z = theim.nplanes * theim.vz / voxel.z;
      if ( BAL_InitAllocImage( &resim, (char*)NULL, dim.x, dim.y, dim.z, 1, theim.type ) != 1 ) {
	BAL_FreeTransformation( &theTrsf );
	BAL_FreeImage( &theim );
	if ( isVerbose )
	  fprintf( stderr, "%s: unable to initialize result image (from parameters, 3D case)\n", program );
	return -1;
      }
    }
    else {
      if ( BAL_InitAllocImage( &resim, (char*)NULL, dim.x, dim.y, 1, 1, theim.type ) != 1 ) {
	BAL_FreeTransformation( &theTrsf );
	BAL_FreeImage( &theim );
	if ( isVerbose )
	  fprintf( stderr, "%s: unable to initialize result image (from parameters, 2D case)\n", program );
	return -1;
      }
    }
    if ( voxel.x > 0.0 ) resim.vx = voxel.x;
    if ( voxel.y > 0.0 ) resim.vy = voxel.y;
    if ( voxel.z > 0.0 ) resim.vz = voxel.z;
  }

  /* initialisation with transformation
   */
  else if ( BAL_IsTransformationVectorField( &theTrsf ) == 1 ) {
    if ( BAL_InitAllocImage( &resim, (char*)NULL, 
			     theTrsf.vx.ncols, theTrsf.vx.nrows, theTrsf.vx.nplanes, 
			     theim.vdim, theim.type ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      BAL_FreeImage( &theim );
      if ( isVerbose )
	fprintf( stderr, "%s: unable to initialize result image (from vector field)\n", program );
      return -1;
    }
    resim.vx = theTrsf.vx.vx;
    resim.vy = theTrsf.vx.vy;
    resim.vz = theTrsf.vx.vz;
  }

  /* initialisation with input image
   */
  else {
    if ( BAL_InitAllocImage( &resim, (char*)NULL, 
			     theim.ncols, theim.nrows, theim.nplanes, 
			     theim.vdim, theim.type ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      BAL_FreeImage( &theim );
      if ( isVerbose )
	fprintf( stderr, "%s: unable to initialize result image (from input image)\n", program );
      return -1;
    }
    resim.vx = theim.vx;
    resim.vy = theim.vy;
    resim.vz = theim.vz;
  }


  /***************************************************
   *
   * 
   *
   ***************************************************/

  /* resampling transformation
     - the given transformation, if any
     - calculate a image-to-image transformation
  */

  if ( real_transformation_name == NULL && voxel_transformation_name == NULL ) {
    if ( BAL_AllocTransformation( &theTrsf, AFFINE_3D, (bal_image *)NULL ) != 1 ) {
      BAL_FreeImage( &resim );
      BAL_FreeImage( &theim );
      if ( isVerbose )
	fprintf( stderr, "%s: unable to initialize transformation\n", program );
      return -1;
    }
    
    if ( resize && (dim.x > 0 && dim.y > 0) ) {
      resim.vx = theim.vx * (float)theim.ncols / (float)resim.ncols;
      resim.vy = theim.vy * (float)theim.nrows / (float)resim.nrows;
      resim.vz = theim.vz * (float)theim.nplanes / (float)resim.nplanes;
    }

    if ( BAL_ComputeImageToImageTransformation( &resim, &theim, &theTrsf ) != 1 ) {
      BAL_FreeImage( &resim );
      BAL_FreeImage( &theim );
      if ( isVerbose )
	fprintf( stderr, "%s: unable to compute transformation\n", program );
      return -1;
    }
  }


  
  /***************************************************
   *
   * 
   *
   ***************************************************/

  if ( isDebug )
    BAL_PrintTransformation( stderr, &theTrsf, "resampling transformation" ) ;

  if ( BAL_ResampleImage( &theim, &resim, &theTrsf, interpolation ) != 1 ) {
    BAL_FreeImage( &resim );
    BAL_FreeTransformation( &theTrsf );
    BAL_FreeImage( &theim );
    if ( isVerbose )
      fprintf( stderr, "%s: unable to compute resampling\n", program );
    return -1;
  }


  if ( result_real_transformation_name != (char*)NULL ) {
    if ( BAL_WriteTransformation( &theTrsf, result_real_transformation_name ) != 1 ) {
      BAL_FreeImage( &resim );
      BAL_FreeTransformation( &theTrsf );
      BAL_FreeImage( &theim );
      if ( isVerbose )
	fprintf( stderr, "%s: unable to write real transformation\n", program );
      return -1;
    }
  }

  if ( result_voxel_transformation_name != (char*)NULL ) {
    if ( BAL_ChangeTransformationToVoxelUnit( &theim, &resim, &theTrsf, &theTrsf ) != 1 ) {
      BAL_FreeImage( &resim );
      BAL_FreeTransformation( &theTrsf );
      BAL_FreeImage( &theim );
      if ( isVerbose )
	fprintf( stderr, "%s: unable to convert transformation\n", program );
      return -1;
    }
    if ( BAL_WriteTransformation( &theTrsf, result_voxel_transformation_name ) != 1 ) {
      BAL_FreeImage( &resim );
      BAL_FreeTransformation( &theTrsf );
      BAL_FreeImage( &theim );
      if ( isVerbose )
	fprintf( stderr, "%s: unable to write voxel transformation\n", program );
      return -1;
    }
  }
 
  BAL_FreeTransformation( &theTrsf );
  BAL_FreeImage( &theim );

  if ( type == resim.type || type == TYPE_UNKNOWN ) {
    
    if ( BAL_WriteImage( &resim, resim_name ) != 1 ) {
      BAL_FreeImage( &resim );
      if ( isVerbose )
	fprintf( stderr, "%s: unable to write result image '%s'\n", program, resim_name );
      return -1;
    }

  }
  else {

    if ( BAL_InitAllocImage( &auxim, (char*)NULL, 
			     resim.ncols, resim.nrows, resim.nplanes, 
			     resim.vdim, type ) != 1 ) {
      BAL_FreeImage( &resim );
      if ( isVerbose )
	fprintf( stderr, "%s: unable to allocate auxiliary result image\n", program );
      return -1;
    }

    if ( BAL_CopyImage( &resim, &auxim ) != 1 ) {
      BAL_FreeImage( &auxim );
      BAL_FreeImage( &resim );
      if ( isVerbose )
	fprintf( stderr, "%s: unable to convert result image\n", program );
      return -1;
    }
    
    if ( BAL_WriteImage( &auxim, resim_name ) != 1 ) {
      BAL_FreeImage( &auxim );
      BAL_FreeImage( &resim );
      if ( isVerbose )
	fprintf( stderr, "%s: unable to write result image '%s'\n", program, resim_name );
      return -1;
    }

    BAL_FreeImage( &auxim );

  }

  BAL_FreeImage( &resim );


  return 0;
}
