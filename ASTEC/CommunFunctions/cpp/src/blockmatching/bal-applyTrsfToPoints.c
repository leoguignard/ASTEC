/*************************************************************************
 * applyTrsf.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2012, all rights reserved
 *
 * AUTHOR:
 * Etienne Delclaux (etienne.delclaux@inria.fr) 
 * From Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Tue Sep 17 12:39:00 CET 2013
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */
#include <bal-applyTrsfToPoints.h>

static char* program = "applyTrsfToPoints";

int applyTrsfToPoints(
       char *thefile_name,
       char *resfile_name,
       char *real_transformation_name,
       int verbose,
       int debug)
{
  bal_transformation theTrsf; 
  FILE *theF;
  FILE *resF;

  bal_doublePoint thePt;
  bal_doublePoint resPt;
  double x, y, z;
  int i, r;


  BAL_InitTransformation( &theTrsf );

  /* reading transformation, if any 
   */
  if ( real_transformation_name != NULL ) {
    if ( BAL_ReadTransformation( &theTrsf, real_transformation_name ) != 1 ) {
      if ( verbose )
	fprintf( stderr, "%s: unable to read 'real' transformation '%s'\n", program, real_transformation_name );
      return -1;
    }
  }
  else {
    if ( verbose )
      fprintf( stderr, "%s: no transformation \n", program );
    return -1;
  }



  theF = fopen( thefile_name, "r" );
  if ( theF == NULL ) {
    BAL_FreeTransformation( &theTrsf );
    if ( verbose )
      fprintf( stderr, "%s: unable to open '%s' for reading\n", program, thefile_name );
    return -1;
  }
  
  resF = fopen( resfile_name, "w" );
  if ( resF == NULL ) {
    fclose( theF );
    BAL_FreeTransformation( &theTrsf );
    if ( verbose )
      fprintf( stderr, "%s: unable to open '%s' for writing\n", program, resfile_name );
    return -1;
  }



  i = 0;
  while ( (r=fscanf( theF, "%lf %lf %lf", &x, &y, &z )) == 3 || r == 2 ) {
    if ( r == 2 ) z = 0.0;

    if ( debug ) 
      fprintf( stderr, "#%04d: read %lf %lf %lf\n", i, x, y, z );

    thePt.x = x;
    thePt.y = y;
    thePt.z = z;
    if ( BAL_TransformPoint( &thePt, &resPt, &theTrsf ) != 1 ) {
      fclose( resF );
      fclose( theF );
      BAL_FreeTransformation( &theTrsf );
      if ( verbose )
	fprintf( stderr, "%s: unable to transform point #%d\n", program, i );
      return -1;
    }

    if ( debug ) 
      fprintf( stderr, "\t write %lf %lf %lf\n", resPt.x, resPt.y, resPt.z );

    fprintf( resF, "%f %f %f\n", resPt.x, resPt.y, resPt.z );
    i++;
  }

  fclose( resF );
  fclose( theF );
  BAL_FreeTransformation( &theTrsf );


  return 0;
}
