/*************************************************************************
 * invTrsf.c -
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

#include "bal-invTrsf.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <bal-stddef.h>
#include <bal-transformation.h>
#include <bal-transformation-tools.h>

static char* program = "invTrsf";

int invTrsf(
	    char* thetrsf_name,
	    char* restrsf_name,
	    char* template_image_name,
	    bal_integerPoint dim,
	    bal_doublePoint voxel,
	    int isVerbose
	    )
{
  bal_transformation theTrsf;
  bal_transformation resTrsf;
  bal_image template;

  /***************************************************
   *
   * 
   *
   ***************************************************/
  BAL_InitTransformation( &theTrsf );
  BAL_InitTransformation( &resTrsf );



  if ( thetrsf_name == NULL ) {
    if ( isVerbose )
      fprintf( stderr, "%s: no input transformation\n", program );
    return -1;
  }


  if ( BAL_ReadTransformation( &theTrsf, thetrsf_name ) != 1 ) {
    if ( isVerbose )
      fprintf( stderr, "%s: unable to read '%s'\n", program, thetrsf_name );
    return -1;
  }



  switch ( theTrsf.type ) {
    
  default :

    if ( isVerbose )
      BAL_PrintTransformation( stderr, &theTrsf, "read transformation" );
    fprintf( stderr, "%s: such transformation type not handled yet\n", program );
    BAL_FreeTransformation( &theTrsf );
    return -1;
    
  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :

    if ( BAL_AllocTransformation( &resTrsf, theTrsf.type, (bal_image *)NULL ) != 1 ) {
      if ( isVerbose )
	fprintf( stderr, "%s: unable to allocate result transformation (linear case)\n", program );
      BAL_FreeTransformation( &theTrsf );
      return -1;
    }
    break;

  case VECTORFIELD_2D :
  case VECTORFIELD_3D :

    /* initializing result image
       - with reference image, if any
    */
    if ( template_image_name != NULL ) {
      if ( BAL_ReadImage( &template, template_image_name, 1 ) != 1 ) {
	if ( isVerbose )
	  fprintf( stderr, "%s: unable to read '%s'\n", program, template_image_name );
	BAL_FreeTransformation( &theTrsf );
	return -1;
      }
    }

    /* initializing result image
       - with parameters, if any
    */
    else if ( dim.x > 0 && dim.y > 0 ) {
      if ( dim.z > 0 ) {
	if ( BAL_InitImage( &template, (char*)NULL, dim.x, dim.y, dim.z, 1, UCHAR ) != 1 ) {
	  if ( isVerbose )
	    fprintf( stderr, "%s: unable to initialize auxiliary image\n", program );
	  BAL_FreeTransformation( &theTrsf );
	  return -1;
	}
      }
      else {
	if ( BAL_InitImage( &template, (char*)NULL, dim.x, dim.y, 1, 1, UCHAR ) != 1 ) {
	  if ( isVerbose )
	    fprintf( stderr, "%s: unable to initialize auxiliary image (dimz=1) \n", program );
	  BAL_FreeTransformation( &theTrsf );
	  return -1;
	}
      }
      if ( voxel.x > 0.0 ) template.vx = voxel.x;
      if ( voxel.y > 0.0 ) template.vy = voxel.y;
      if ( voxel.z > 0.0 ) template.vz = voxel.z;
    }
    
    /* initialisation with transformation
     */
    else {
      if ( BAL_InitImage( &template, (char*)NULL, 
			  theTrsf.vx.ncols, theTrsf.vx.nrows, theTrsf.vx.nplanes, 
			  1, UCHAR ) != 1 ) {
	BAL_FreeTransformation( &theTrsf );
	fprintf( stderr, "%s: unable to initialize  auxiliary image (input vector field)\n", program );
	return -1;
      }
      template.vx = theTrsf.vx.vx;
      template.vy = theTrsf.vx.vy;
      template.vz = theTrsf.vx.vz;
    }
    
    if ( theTrsf.type == VECTORFIELD_2D ) {
      template.nplanes = theTrsf.vx.nplanes;
    }
    
    
    if ( BAL_AllocTransformation( &resTrsf, theTrsf.type, &template ) != 1 ) {
      if ( isVerbose )
	fprintf( stderr, "%s: unable to allocate result transformation (vector field)\n", program );
      BAL_FreeImage( &template );
      BAL_FreeTransformation( &theTrsf );
      return -1;
    }

    BAL_FreeImage( &template );

  }



  if ( BAL_InverseTransformation(  &theTrsf, &resTrsf ) != 1 ) {
    BAL_FreeTransformation( &theTrsf );
    BAL_FreeTransformation( &resTrsf );
    if ( isVerbose )
      fprintf( stderr, "%s: unable to invert transformation '%s'\n", program, thetrsf_name );
    return -1;
  }





  BAL_FreeTransformation( &theTrsf );
  
  if ( restrsf_name != NULL ) {
    if ( BAL_WriteTransformation( &resTrsf, restrsf_name ) != 1 ) {
      if ( isVerbose )
	fprintf( stderr, "%s: unable to write '%s'\n", program, restrsf_name );
      return -1;
    }
  }
  else {
    if ( isVerbose )
      fprintf( stderr, "%s: no output transformation name\n", program );
    return -1;
  }


  BAL_FreeTransformation( &resTrsf );


  return 0;
}
