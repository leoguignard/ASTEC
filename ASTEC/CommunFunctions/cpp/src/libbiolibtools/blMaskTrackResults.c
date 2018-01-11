/*************************************************************************
 * blMaskTrackResults.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2014, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Lun 17 fév 2014 14:16:59 CET
 *
 * ADDITIONS, CHANGES
 *
 */

#include <stdlib.h>
#include <stdio.h>

#include <sys/time.h> /* gettimeofday() */
#include <time.h> /* clock() */
#include <string.h>

static int _verbose_ = 1;
static int _debug_ = 0;


#include <bal-image.h>
#include <bal-biolib-tools.h>






typedef struct local_parameter {

  /* file names
   */
  char *thetrack_name;
  char *restrack_name;

  char *themask_name;

} local_parameter;





/*------- Definition des fonctions statiques ----------*/
static void _ErrorParse( char *str, int flag );
static void _Parse( int argc, char *argv[], local_parameter *p );
static void _InitParam( local_parameter *par );




static char *program = NULL;

static char *usage = "%s %s [-mask |-m %s]\n\
 [-v] [-help]";

static char *detail = "\
-mask|-m %s # mask to be applied\n\
 -v : mode verbose\n\
\n";






static int _process( char *thetrack_name,
		     char *restrack_name,
		     char *themask_name );





int main( int argc, char *argv[] )
{
  local_parameter p;

  /***************************************************
   *
   * parsing parameters
   *
   ***************************************************/
  program = argv[0];

  /* no arguments
   */
  if ( argc == 1 ) _ErrorParse( NULL, 0 );


  /* parsing parameters 
   */
  _InitParam( &p );
  _Parse( argc, argv, &p );
   


  /****************************************
   *
   * single file case
   *
   ****************************************/
  if ( p.thetrack_name != (char*)NULL
       && p.themask_name != (char*)NULL
       && p.restrack_name != (char*)NULL ) {

    if ( _process( p.thetrack_name, p.restrack_name,
		   p.themask_name ) != 0 ) {
      _ErrorParse( "unable to process single file case ...\n", 0);
    }
    
  }


  /****************************************
   *
   * undefined case
   *
   ****************************************/
  else {
    _ErrorParse( "I don't know what to do ...\n", 0);
  }



  return( 1 );
}







static int _process( char *thetrack_name,
		     char *restrack_name,
		     char *themask_name ) 
{
  char *proc = "_process";
  bal_blTrackList thelist, reslist;
  bal_image theMask;
  int n, i;
  int x, y, z;
  int isincluded;

  if ( 0 ) {
    fprintf( stderr, "%s: processing '%s' into '%s' with '%s'\n",
	     proc, thetrack_name,
	     restrack_name,
	     themask_name );
  }

  /* track lists
   */
  BAL_InitBlTrackList( &thelist );
  if ( BAL_ReadBlTrackList( &thelist, thetrack_name ) != 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to read track file '%s'\n", proc, thetrack_name );
    return( -1 );
  }
  
  BAL_InitBlTrackList( &reslist );
  if ( BAL_AllocBlTrackList( &reslist, thelist.n ) != 0 ) {
    BAL_FreeBlTrackList( &thelist );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate result track\n", proc );
    return( -1 );
  }
  reslist.n = thelist.n;


  /* reading mask 
   */
  BAL_InitImage( &theMask, NULL, 0, 0, 0, 0, UCHAR );
  if ( BAL_ReadImage( &theMask, themask_name, 0 ) != 1 ) {
    _ErrorParse( "error when reading mask image\n", 0);
  }


  /* track processing
   */
  switch ( theMask.type ) {
  default :
    BAL_FreeImage( &theMask );
    BAL_FreeBlTrackList( &reslist );
    BAL_FreeBlTrackList( &thelist );
    if ( _verbose_ )
      fprintf( stderr, "%s: mask image type not handled yet\n", 
	       proc );
    return( -1 );
    break;
  case UCHAR :
    {
      u8 ***theBuf = (u8***)theMask.array;
      for ( n=0, reslist.n=0; n<thelist.n; n++ ) {
	for ( isincluded=1, i=0; i<thelist.data[n].detectionList.n && isincluded==1; i++ ) {
	  x = (int)(thelist.data[n].detectionList.data[i].voxelcenter.x + 0.5);
	  y = (int)(thelist.data[n].detectionList.data[i].voxelcenter.y + 0.5);
	  z = (int)(thelist.data[n].detectionList.data[i].voxelcenter.z + 0.5);
	  if ( x < 0 || x >= theMask.ncols ) isincluded = 0;
	  if ( y < 0 || y >= theMask.nrows ) isincluded = 0;
	  if ( z < 0 || z >= theMask.nplanes ) isincluded = 0;
	  if ( theBuf[z][y][x] == 0 ) isincluded = 0;
	}
	if ( isincluded == 1 ) {
	  if ( BAL_CopyBlTrack( &(thelist.data[n]),  &(reslist.data[reslist.n]) ) != 1 ) {
	    BAL_FreeBlTrackList( &reslist );
	    BAL_FreeBlTrackList( &thelist );
	    if ( _verbose_ )
	      fprintf( stderr, "%s: unable to copy track #%d\n", proc, n );
	    return( -1 );
	  }
	  reslist.n ++;
	}
      }
    }
    break;
  }

  BAL_FreeImage( &theMask );
  BAL_FreeBlTrackList( &thelist );
  


  /* result
   */
  
  if ( BAL_WriteBlTrackList( &reslist, restrack_name ) != 0 ) {
    BAL_FreeBlTrackList( &reslist );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to write track file '%s'\n", proc, restrack_name );
    return( -1 );
  }
  
  BAL_FreeBlTrackList( &reslist );

  return( 0 );
}










static void _Parse( int argc, char *argv[], local_parameter *p )
{
  int i;
  int inputisread = 0;
  int outputisread = 0;

  program = argv[0];
	
  for ( i=1; i<argc; i++ ) {
  
    if ( argv[i][0] == '-' ) {

      /* mask related file names
       */

      if ( strcmp ( argv[i], "-mask") == 0
	   || (strcmp ( argv[i], "-m") == 0 && argv[i][2] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -mask", 0 );
	p->themask_name = argv[i];
      }
      
      /* ...
       */


      
      
      else if ( (strcmp ( argv[i], "-verbose") == 0 && argv[i][8] == '\0')
		|| (strcmp ( argv[i], "-v") == 0 && argv[i][2] == '\0') ) {
	if ( _verbose_ <= 0 ) _verbose_ = 1;
	else _verbose_ ++;
      }
      
      else if ( strcmp ( argv[i], "--help" ) == 0 
		|| ( strcmp ( argv[i], "-help" ) == 0 && argv[i][5] == '\0' )
		|| ( strcmp ( argv[i], "--h" ) == 0 && argv[i][3] == '\0' )
		|| ( strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0' ) ) {
	_ErrorParse( NULL, 1 );
      }
      
       else if ( (strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0') ) {
	_debug_ = 1;
      }

      /* unknown option
       */
      else {
	fprintf(stderr,"unknown option: '%s'\n",argv[i]);
      }
    }
    
    /*--- saisie des noms d'images ---*/
    else if ( argv[i][0] != 0 ) {
      if ( inputisread == 0 ) {
	p->thetrack_name = argv[i];
	inputisread = 1;
      }
      else if ( outputisread == 0 ) {
	p->restrack_name = argv[i];
	outputisread = 1;
      }
      else 
	fprintf(stderr,"too many file names: '%s'\n",argv[i]);
    }

  }
  
}





static void _ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  if ( str != (char*)NULL )
    (void)fprintf(stderr,"Erreur : %s",str);
  exit(0);
}





static void _InitParam( local_parameter *p )
{
  p->thetrack_name = (char*)NULL;
  p->restrack_name = (char*)NULL;

  p->themask_name = (char*)NULL;
}
