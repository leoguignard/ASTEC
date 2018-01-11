/*************************************************************************
 * bal-createRandomTrsf.c -
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


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <bal-transformation-tools.h>



#include <bal-createRandomTrsf.h>



static char *program = "createRandomTrsf";
static int _verbose_ = 1;

int createRandomTrsf(
        char *restrsf_name,
        char *template_name,
        bal_doublePoint fixedpoint,
        enumTypeTransfo transformation_type,
        int print)
{
  bal_transformation theTrsf;
  long int seedRandom = time(0);
  bal_image templateImage;
  srandom( seedRandom );

  /***************************************************
   *
   * 
   *
   ***************************************************/
  BAL_InitTransformation( &theTrsf );
  switch ( transformation_type ) {

  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such type not handled yet for input transformation\n", program );
    return -1;
    
  case TRANSLATION_2D :
    if ( BAL_Random2DTranslationMatrix( &theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 2D translation matrix\n", program );
      return -1;
    }
    break;

  case TRANSLATION_3D :
    if ( BAL_Random3DTranslationMatrix( &theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 3D translation matrix\n", program );
      return -1;
    }
    break;

  case TRANSLATION_SCALING_2D :
    if ( BAL_Random2DTranslationScalingMatrix( &theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 2D translation and scaling matrix\n", program );
      return -1;
    }
    break;

  case TRANSLATION_SCALING_3D :
    if ( BAL_Random3DTranslationScalingMatrix( &theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 3D translation and scaling matrix\n", program );
      return -1;
    }
    break;

  case RIGID_2D :
    if ( BAL_Random2DRigidMatrix( &theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 2D rigid matrix\n", program );
      return -1;
    }
    break;

  case RIGID_3D :
    if ( BAL_Random3DRigidMatrix( &theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 3D rigid matrix\n", program );
      return -1;
    }
    break;

  case SIMILITUDE_2D :
    if ( BAL_Random2DSimilitudeMatrix( &theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 2D similitude matrix\n", program );
      return -1;
    }
    break;

  case SIMILITUDE_3D :
    if ( BAL_Random3DSimilitudeMatrix( &theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 3D similitude matrix\n", program );
      return -1;
    }
    break;

  case AFFINE_2D :
    if ( BAL_Random2DAffineMatrix( &theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 2D affine matrix\n", program );
      return -1;
    }
    break;

  case AFFINE_3D :
    if ( BAL_Random3DAffineMatrix( &theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 3D affine matrix\n", program );
      return -1;
    }
    break;
  }

  /* transformation centering
   */
  if ( template_name != NULL ) {
    if ( BAL_ReadImage( &templateImage, template_name, 0 ) != 1 ) {
      fprintf( stderr, "%s: can not read '%s'\n", program, template_name );
      return -1;
    }
    fixedpoint.x = ((templateImage.ncols-1) * templateImage.vx) / 2.0;
    fixedpoint.y = ((templateImage.nrows-1) * templateImage.vy) / 2.0;
    fixedpoint.z = ((templateImage.nplanes-1) * templateImage.vz) / 2.0;
    BAL_FreeImage( &templateImage );
  }

  if ( fixedpoint.x > -9000 && fixedpoint.y > -9000 && fixedpoint.z > -9000 ) {
    switch ( transformation_type ) {
    case VECTORFIELD_2D :
    case VECTORFIELD_3D :
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such type not handled yet for input transformation\n", program );
      return -1;
      
    case TRANSLATION_2D :
    case TRANSLATION_SCALING_2D :
    case RIGID_2D :
    case SIMILITUDE_2D :
    case AFFINE_2D :
      theTrsf.mat.m[3] = fixedpoint.x - theTrsf.mat.m[0] * fixedpoint.x - theTrsf.mat.m[1] * fixedpoint.y;
      theTrsf.mat.m[7] = fixedpoint.y - theTrsf.mat.m[4] * fixedpoint.x - theTrsf.mat.m[5] * fixedpoint.y;
      break;
    case TRANSLATION_3D :
    case TRANSLATION_SCALING_3D :
    case RIGID_3D :
    case SIMILITUDE_3D :
    case AFFINE_3D :
      theTrsf.mat.m[ 3] = fixedpoint.x - theTrsf.mat.m[0] * fixedpoint.x - theTrsf.mat.m[1] * fixedpoint.y - theTrsf.mat.m[ 2] * fixedpoint.z;
      theTrsf.mat.m[ 7] = fixedpoint.y - theTrsf.mat.m[4] * fixedpoint.x - theTrsf.mat.m[5] * fixedpoint.y - theTrsf.mat.m[ 6] * fixedpoint.z;
      theTrsf.mat.m[11] = fixedpoint.z - theTrsf.mat.m[8] * fixedpoint.x - theTrsf.mat.m[9] * fixedpoint.y - theTrsf.mat.m[10] * fixedpoint.z;
      break;
    }
  }

  if ( print )
    BAL_PrintTransformation( stderr, &theTrsf, "generated transformation" );

  if ( restrsf_name != NULL ) {
    if ( BAL_WriteTransformation( &theTrsf, restrsf_name ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write '%s'\n", program, restrsf_name );
      BAL_FreeTransformation( &theTrsf );
      return -1;
    }
  }

  BAL_FreeTransformation( &theTrsf );
  return 1;
}
