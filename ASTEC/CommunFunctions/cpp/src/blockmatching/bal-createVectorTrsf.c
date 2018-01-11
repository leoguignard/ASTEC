/*************************************************************************
 * createRandomTrsf.c -
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

#include "bal-createVectorTrsf.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <bal-transformation-tools.h>

static char *program = "createVectorTrsf";

int createVectorTrsf(
		     char *restrsf_name,
		     char *template_image_name,
		     bal_integerPoint dim,
		     bal_doublePoint voxel,
		     enumTypeTransfo transformation_type,
		     enumVectorType vector_type,
		     int isDebug,
		     int isVerbose
		     )
{
  bal_transformation theTrsf;
  bal_image template;

  /***************************************************
   *
   * 
   *
   ***************************************************/
  BAL_InitTransformation( &theTrsf );

  switch ( transformation_type ) {

  default :
    if ( isVerbose )
      fprintf( stderr, "%s: such type not handled yet for input transformation\n", program );
    return( -1 );

  case VECTORFIELD_2D :
  case VECTORFIELD_3D :

    /* initializing result image
       - with reference image, if any
    */
    if ( template_image_name != NULL ) {
      if ( BAL_ReadImage( &template, template_image_name, 1 ) != 1 ) {
	if ( isVerbose )
	  fprintf( stderr, "%s: unable to read '%s'\n", program, template_image_name );
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
	  return -1;
	}
      }
      else {
	if ( BAL_InitImage( &template, (char*)NULL, dim.x, dim.y, 1, 1, UCHAR ) != 1 ) {
	  if ( isVerbose )
	    fprintf( stderr, "%s: unable to initialize auxiliary image (dimz=1) \n", program );
	  return -1;
	}
      }
      if ( voxel.x > 0.0 ) template.vx = voxel.x;
      if ( voxel.y > 0.0 ) template.vy = voxel.y;
      if ( voxel.z > 0.0 ) template.vz = voxel.z;
    }
    else {
      if ( isVerbose )
	fprintf( stderr, "%s: negative dimensions, unable to initialize auxiliary image\n", program );
      return -1;
    }
    
    if ( BAL_AllocTransformation( &theTrsf, transformation_type, &template ) != 1 ) {
      if ( isVerbose )
	fprintf( stderr, "%s: unable to allocate result transformation\n", program );
      BAL_FreeImage( &template );
      return -1;
    }

    BAL_FreeImage( &template );

  }



  switch ( vector_type ) {			

  default :
    BAL_FreeTransformation( &theTrsf );
    if ( isVerbose )
      fprintf( stderr, "%s: such type not handled yet for vector field\n", program );
    return( -1 );

  case SINUS_2D :
    if (  BAL_Sinusoid2DVectorField( &theTrsf ) != 1 ) {
      if ( isVerbose )
	fprintf( stderr, "%s: unable to generate 2D sinus transformation\n", program );
      return -1;
    }
    break;

  case SINUS_3D :
    if (  BAL_Sinusoid3DVectorField( &theTrsf ) != 1 ) {
      if ( isVerbose )
	fprintf( stderr, "%s: unable to generate 3D sinus transformation\n", program );
      return -1;
    }
    break;
    


  }




  if ( isDebug )
    BAL_PrintTransformation( stderr, &theTrsf, "generated transformation" );

  if ( restrsf_name != NULL ) {
    if ( BAL_WriteTransformation( &theTrsf, restrsf_name ) != 1 ) {
      if ( isVerbose )
	fprintf( stderr, "%s: unable to write '%s'\n", program, restrsf_name );
      BAL_FreeTransformation( &theTrsf );
      return -1;
    }
  }
    

  BAL_FreeTransformation( &theTrsf );


  return 0;
}
