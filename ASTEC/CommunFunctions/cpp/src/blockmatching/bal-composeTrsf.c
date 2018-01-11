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
#include "bal-composeTrsf.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <bal-stddef.h>
#include <bal-image.h>
#include <bal-image-tools.h>
#include <bal-transformation.h>
#include <bal-transformation-tools.h>

static char *program = "composeTrsf";

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
	int isDebug,
	int isVerbose
        )
{
  /* initialize parameters
   */
  int i;
  bal_transformation theTrsf;
  bal_transformation *resTrsf;
  bal_transformation tmp1Trsf;
  bal_transformation tmp2Trsf;
  bal_image template;

  BAL_InitTransformation( &theTrsf );
  BAL_InitTransformation( &tmp1Trsf );
  BAL_InitTransformation( &tmp2Trsf );

  /***************************************************
   *
   * 
   *
   ***************************************************/
  /* find the first transformation
   */
  for ( i=argc-1; i>=0 && is_a_trsf[i] == 0; i-- )
    ;

  if ( i <= 0 ) {
    fprintf( stderr, "no transformations to be composed ? (use -trsfs)");
  }

  if ( isDebug ) {
    fprintf( stderr, " first transformation          = argc %d = '%s'\n", i, argv[i] );
  }

  if ( BAL_ReadTransformation( &tmp1Trsf, argv[i] ) != 1 ) {
    if ( isVerbose )
      fprintf( stderr, "%s: unable to read '%s'\n", program, argv[i] );
   return -1;
  }
  resTrsf = &tmp1Trsf;

  

  /* loop over transformations
   */
  for ( i-- ; i>=0; i-- ) {
    if ( is_a_trsf[i] == 0 ) continue;
     if ( isDebug ) {
       fprintf( stderr, " transformation to be composed = argc %d = '%s'\n", i, argv[i] );
     }
     
     if ( BAL_ReadTransformation( &theTrsf, argv[i] ) != 1 ) {
       BAL_FreeTransformation(  resTrsf );
       if ( isVerbose )
	 fprintf( stderr, "%s: unable to read '%s'\n", program, argv[i] );
      return -1;
     }



     /* check whether the current transformation is a matrix
	and the next one is a vector field
	if yes, change the current into a vector field
     */

     if ( BAL_IsTransformationLinear( resTrsf ) == 1 ) {

       if ( BAL_IsTransformationVectorField( &theTrsf ) == 1 ) {

	 /* initializing result transformation
	    - with reference image, if any
	    - with parameters, if any
	    - with transformation
	 */

	 /* initialisation
	    - with reference image, if any
	  */
	 if ( template_image_name != NULL ) {
	   if ( BAL_ReadImage( &template, template_image_name, 0 ) != 1 ) {
	     BAL_FreeTransformation(  resTrsf );
	     fprintf( stderr, "%s: can not read reference image '%s'\n", program, template_image_name );
	    return -1;
	   }
	 }

	 /* initialisation
	    - with parameters, if any
	 */
	 else if ( dim.x > 0 && dim.y > 0 ) {
	   if ( dim.z > 0 ) {
	     if ( BAL_InitImage( &template, (char*)NULL, dim.x, dim.y, dim.z, 1, UCHAR ) != 1 ) {
	       BAL_FreeTransformation(  resTrsf );
	       fprintf( stderr, "%s: unable to initialize result transformation (parameters, 3D case)\n", program );
	      return -1;
	     }
	   }
	   else {
	     if ( BAL_InitImage( &template, (char*)NULL, dim.x, dim.y, 1, 1, UCHAR ) != 1 ) {
	       BAL_FreeTransformation(  resTrsf );
	       fprintf( stderr, "%s: unable to initialize result transformation (parameters, 2D case)\n", program );
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
	     BAL_FreeTransformation(  resTrsf );
	     fprintf( stderr, "%s: unable to initialize result transformation (vector field)\n", program );
	    return -1;
	   }
	   template.vx = theTrsf.vx.vx;
	   template.vy = theTrsf.vx.vy;
	   template.vz = theTrsf.vx.vz;
	 }
	 

	 /* allocation of a new auxiliary transformation
	  */
	 if ( BAL_AllocTransformation( &tmp2Trsf, theTrsf.type, &template ) != 1 ) {
	   BAL_FreeTransformation(  resTrsf );
	   fprintf( stderr, "%s: unable to allocate new auxiliary transformation\n", program );
	  return -1;
	 }

	 if ( BAL_CopyTransformation( resTrsf, &tmp2Trsf ) != 1 ) {
	   BAL_FreeTransformation( &tmp2Trsf );
	   BAL_FreeTransformation( resTrsf );
	   fprintf( stderr, "%s: unable to copy transformation\n", program );
	  return -1;
	 }

	 BAL_FreeTransformation( resTrsf );
	 resTrsf =  &tmp2Trsf;

       }
     }




     /* compose transformations
	int BAL_TransformationComposition( bal_transformation *t_res,
	                               bal_transformation *t1, 
	                               bal_transformation *t2 ) 
	with t_res = t1 o t2 
     */
     if ( BAL_TransformationComposition( resTrsf, &theTrsf, resTrsf ) != 1 ) {
       BAL_FreeTransformation(  &theTrsf );
       BAL_FreeTransformation(  resTrsf );
       if ( isVerbose )
	 fprintf( stderr, "%s: unable to compose '%s' with intermediary result\n", program, argv[i] );
      return -1;
     }

     /* free the additional transformation
	for next reading
     */
     BAL_FreeTransformation(  &theTrsf );

  }

  if ( BAL_WriteTransformation( resTrsf, restrsf_name ) != 1 ) {
    if ( isVerbose )
      fprintf( stderr, "%s: unable to write '%s'\n", program, restrsf_name );
    BAL_FreeTransformation( resTrsf );
   return -1;
  }
  
  BAL_FreeTransformation( resTrsf );

  return 0;
}
