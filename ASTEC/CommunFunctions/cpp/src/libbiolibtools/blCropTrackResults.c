/*************************************************************************
 * blCropTrackResults.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2014, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Lun 17 fév 2014 12:06:41 CET
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
static int _debug_ = 1;

#include <bal-image.h>
#include <bal-transformation.h>
#include <bal-biolib-tools.h>






typedef struct local_parameter {

  /* file names
   */
  char *thetrack_name;
  char *restrack_name;

  /* crop parameters
   */
  char *template_image_name; /* template */  

  bal_integerPoint origin;
  
  bal_integerPoint dim;
  
  bal_integerPoint slice;

} local_parameter;





/*------- Definition des fonctions statiques ----------*/
static void _ErrorParse( char *str, int flag );
static void _Parse( int argc, char *argv[], local_parameter *p );
static void _InitParam( local_parameter *par );




static char *program = NULL;

static char *usage = "%s %s\n\
 [-origin %d %d [%d] | [-ix %d] [-iy %d] [-iz %d]]\n\
 [-dim %d %d [%d] | [-x %d] [-y %d] [-z %d]] [-template %s]\n\
 [-xy %d | -yz %d | -xz %d]\n\
 [-v] [-help]";

static char *detail = "\
-origin %d %d [%d]   # origin of the output image in the input image\n\
                       default is (0,0,0)\n\
-dim %d %d [%d]      # output image dimensions\n\
-template %s         # template image for the dimensions\n\
                       of the output image\n\
-xy %d   # extract the XY slice #%d\n\
         # this is equivalent to '-origin 0 0 %d -dim dimx dimy 1'\n\
-xz %d   # extract the XZ slice #%d\n\
         # the output image has sizes (dimx, dimz, 1)\n\
-yz %d   # extract the YZ slice #%d\n\
         # the output image has sizes (dimy, dimz, 1)\n\
 -v : mode verbose\n\
\n";






static int _process( char *thetrack_name,
		     char *restrack_name,
		     int *leftCorner,
		     int *rightCorner );





int main( int argc, char *argv[] )
{
  local_parameter p;

  bal_image template;
  int theLeftCorner[3] = {-1, -1, -1};
  int theRightCorner[3] = {-1, -1, -1};


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
  

  
  /* initializing result image
     - with reference image, if any
     - with parameters, if any
     - with transformation, if vector field
     - with input image

     - selected coordinates verify: left <= x < right
  */

  /* initialisation with a reference image 
   */
  if ( p.template_image_name != NULL ) {

    if ( BAL_ReadImage( &template, p.template_image_name, 0 ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: can not read template image '%s'\n", program, p.template_image_name );
      return -1;
    }

    theLeftCorner[0] = p.origin.x;
    theLeftCorner[1] = p.origin.y;
    theLeftCorner[2] = p.origin.z;

    theRightCorner[0] = p.origin.x + template.ncols;
    theRightCorner[1] = p.origin.y + template.nrows;
    theRightCorner[2] = p.origin.z + template.nplanes;

    BAL_FreeImage( &template );

    if ( theLeftCorner[0] < 0 ) theLeftCorner[0] = 0;
    if ( theLeftCorner[1] < 0 ) theLeftCorner[1] = 0;
    if ( theLeftCorner[2] < 0 ) theLeftCorner[2] = 0;
    
  }



  /* initialisation with dimension parameters
   */
  else if ( p.dim.x > 0 && p.dim.y > 0 ) {
    
    theLeftCorner[0] = p.origin.x;
    theLeftCorner[1] = p.origin.y;

    theRightCorner[0] = p.origin.x + p.dim.x;
    theRightCorner[1] = p.origin.y + p.dim.y;

    if ( theLeftCorner[0] < 0 ) theLeftCorner[0] = 0;
    if ( theLeftCorner[1] < 0 ) theLeftCorner[1] = 0;
    
    if ( p.dim.z > 0 ) {
      
      theLeftCorner[2] = p.origin.z;
      theRightCorner[2] = p.origin.z + p.dim.z;

      if ( theLeftCorner[2] < 0 ) theLeftCorner[2] = 0;
      
    }
    else {

      theLeftCorner[2] = -1;
      theRightCorner[2] = -1;

    }

  }
  


  /* initialisation with slice information
   */
  else if ( p.slice.z >= 0 ) {

    theLeftCorner[0] = -1;
    theRightCorner[0] = -1;

    theLeftCorner[1] = -1;
    theRightCorner[1] = -1;

    theLeftCorner[2] = p.slice.z;
    theRightCorner[2] = p.slice.z + 1 ;

  }

  else if ( p.slice.y >= 0 ) {

    theLeftCorner[0] = -1;
    theRightCorner[0] = -1;
    
    theLeftCorner[1] = p.slice.y;
    theRightCorner[1] = p.slice.y + 1 ;

    theLeftCorner[2] = -1;
    theRightCorner[2] = -1;

  }

  else if ( p.slice.x >= 0 ) {

    theLeftCorner[0] = p.slice.x;
    theRightCorner[1] = p.slice.x + 1 ;

    theLeftCorner[1] = -1;
    theRightCorner[1] = -1;

    theLeftCorner[2] = -1;
    theRightCorner[2] = -1;

  }

  else {
    if ( _verbose_ )
      fprintf( stderr, "%s: no (valid?) geometrical information?\n", program );
    return -1;
  }

  


  /****************************************
   *
   * single file case
   *
   ****************************************/
  if ( p.thetrack_name != (char*)NULL
       && p.restrack_name != (char*)NULL ) {

    if ( _process( p.thetrack_name, p.restrack_name,
		   theLeftCorner, theRightCorner ) != 0 ) {
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
		     int *leftCorner,
		     int *rightCorner )
{
  char *proc = "_process";
  bal_blTrackList thelist, reslist;
  int i, n;

  if ( 0 ) {
    fprintf( stderr, "%s: processing '%s' into '%s'\n",
	     proc, thetrack_name,
	     restrack_name );
    fprintf( stderr, "   with left=[%d,%d,%d] and right=[%d,%d,%d]\n",
	     leftCorner[0], leftCorner[1], leftCorner[2],
	     rightCorner[0], rightCorner[1], rightCorner[2] );
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



  /* track processing
   */
  for ( n=0, reslist.n=0; n<thelist.n; n++ ) {

    
    BAL_BlTrackBoundingBox( &(thelist.data[n]) );

    if ( leftCorner[0] >= 0 && rightCorner[0] >= 0 ) {
      if ( (int)(thelist.data[n].lcorner.x + 0.5) < leftCorner[0] )
	continue;
      if ( (int)(thelist.data[n].rcorner.x + 0.5) >= rightCorner[0] )
	continue;
    }

    if ( leftCorner[1] >= 0 && rightCorner[1] >= 0 ) {
      if ( (int)(thelist.data[n].lcorner.y + 0.5) < leftCorner[1] )
	continue;
      if ( (int)(thelist.data[n].rcorner.y + 0.5) >= rightCorner[1] )
	continue;
    }

    if ( leftCorner[2] >= 0 && rightCorner[2] >= 0 ) {
      if ( (int)(thelist.data[n].lcorner.z + 0.5) < leftCorner[2] )
	continue;
      if ( (int)(thelist.data[n].rcorner.z + 0.5) >= rightCorner[2] )
	continue;
    }
    
    if ( _debug_ ) {
      fprintf( stderr, "%s: copy input track #%d into output track #%d\n",
	       proc, n, reslist.n );
    }

    if ( BAL_CopyBlTrack( &(thelist.data[n]),  &(reslist.data[reslist.n]) ) != 1 ) {
      BAL_FreeBlTrackList( &reslist );
      BAL_FreeBlTrackList( &thelist );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to copy track #%d\n", proc, n );
      return( -1 );
    }
    
    for ( i=0; i<reslist.data[reslist.n].detectionList.n; i++ ) {
      reslist.data[reslist.n].detectionList.data[i].voxelcenter.x -= leftCorner[0];
      reslist.data[reslist.n].detectionList.data[i].voxelcenter.y -= leftCorner[1];
      reslist.data[reslist.n].detectionList.data[i].voxelcenter.z -= leftCorner[2];
    }

    reslist.n ++;
  }

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
  int i, status;
  int inputisread = 0;
  int outputisread = 0;
  float t;

  program = argv[0];
	
  for ( i=1; i<argc; i++ ) {
  
    if ( argv[i][0] == '-' ) {

      if ( strcmp ( argv[i], "-template") == 0
	   || (strcmp ( argv[i], "-t") == 0 && argv[i][2] == '\0')
	   || (strcmp ( argv[i], "-dims") == 0 && argv[i][5] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -template", 0 );
	p->template_image_name = argv[i];
      }

      /* cropping options
       */
      else if ( (strcmp (argv[i], "-origin" ) == 0 && argv[i][7] == '\0')
		|| (strcmp (argv[i], "-o" ) == 0 && argv[i][2] == '\0') ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -origin %d", 0 );
	status = sscanf( argv[i], "%f", &t );
	if ( status <= 0 ) _ErrorParse( "parsing -origin %d", 0 );
	p->origin.x = (t > 0.0) ? (int)(t+0.5) : (int)(t-0.5);
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -origin %d %d", 0 );
	status = sscanf( argv[i], "%f", &t );
	if ( status <= 0 ) _ErrorParse( "parsing -origin %d %d", 0 );
	p->origin.y = (t > 0.0) ? (int)(t+0.5) : (int)(t-0.5);
	i ++;
	if ( i >= argc) p->origin.z = 0;
	else {
	  status = sscanf( argv[i], "%f", &t );
	  if ( status <= 0 ) {
	    i--;
	    p->origin.z = 0;
	  }
	  else {
	    p->origin.z = (t > 0.0) ? (int)(t+0.5) : (int)(t-0.5);
	  }
	}
      }
      
      else if ( strcmp (argv[i], "-ix" ) == 0 && argv[i][3] == '\0' ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -ix %d", 0 );
	status = sscanf( argv[i], "%f", &t );
	if ( status <= 0 ) _ErrorParse( "parsing -ix %d", 0 );
	p->origin.x = (t > 0.0) ? (int)(t+0.5) : (int)(t-0.5);
      }
      else if ( strcmp (argv[i], "-iy" ) == 0 && argv[i][3] == '\0' ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -iy %d", 0 );
	status = sscanf( argv[i], "%f", &t );
	if ( status <= 0 ) _ErrorParse( "parsing -iy %d", 0 );
	p->origin.y = (t > 0.0) ? (int)(t+0.5) : (int)(t-0.5);
      }
      else if ( strcmp (argv[i], "-iz" ) == 0 && argv[i][3] == '\0' ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -iz %d", 0 );
	status = sscanf( argv[i], "%f", &t );
	if ( status <= 0 ) _ErrorParse( "parsing -iz %d", 0 );
	p->origin.z = (t > 0.0) ? (int)(t+0.5) : (int)(t-0.5);
      }

      else if ( strcmp (argv[i], "-dim" ) == 0 && argv[i][4] == '\0' ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -dim %d", 0 );
	status = sscanf( argv[i], "%d", &(p->dim.x) );
	if ( status <= 0 ) _ErrorParse( "parsing -dim %d", 0 );
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -dim %d %d", 0 );
	status = sscanf( argv[i], "%d", &(p->dim.y) );
	if ( status <= 0 ) _ErrorParse( "parsing -dim %d %d", 0 );
	i ++;
	if ( i >= argc) p->dim.z = 1;
	else {
	  status = sscanf( argv[i], "%d", &(p->dim.z) );
	  if ( status <= 0 ) {
	    i--;
	    p->dim.z = 1;
	  }
	}
      }

      else if ( strcmp (argv[i], "-x" ) == 0 && argv[i][2] == '\0' ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -x %d", 0 );
	status = sscanf( argv[i], "%d", &(p->dim.x) );
	if ( status <= 0 ) _ErrorParse( "parsing -x %d", 0 );
      }
      else if ( strcmp (argv[i], "-y" ) == 0 && argv[i][2] == '\0' ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -y %d", 0 );
	status = sscanf( argv[i], "%d", &(p->dim.y) );
	if ( status <= 0 ) _ErrorParse( "parsing -y %d", 0 );
      }
      else if ( strcmp (argv[i], "-z" ) == 0 && argv[i][2] == '\0' ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -z %d", 0 );
	status = sscanf( argv[i], "%d", &(p->dim.z) );
	if ( status <= 0 ) _ErrorParse( "parsing -z %d", 0 );
      }


      else if ( strcmp (argv[i], "-xy" ) == 0 && argv[i][3] == '\0' ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -xy %d", 0 );
	status = sscanf( argv[i], "%d", &(p->slice.z) );
	if ( status <= 0 ) _ErrorParse( "parsing -xy %d", 0 );
      }
      else if ( strcmp (argv[i], "-yz" ) == 0 && argv[i][3] == '\0' ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -yz %d", 0 );
	status = sscanf( argv[i], "%d", &(p->slice.x) );
	if ( status <= 0 ) _ErrorParse( "parsing -yz %d", 0 );
      }
      else if ( strcmp (argv[i], "-xz" ) == 0 && argv[i][3] == '\0' ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -xz %d", 0 );
	status = sscanf( argv[i], "%d", &(p->slice.y) );
	if ( status <= 0 ) _ErrorParse( "parsing -xz %d", 0 );
      }

      
      
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

  p->template_image_name = (char*)NULL;
  
  p->origin.x = 0;
  p->origin.y = 0;
  p->origin.z = 0;

  p->dim.x = -1;
  p->dim.y = -1;
  p->dim.z = -1;

  p->slice.x = -1;
  p->slice.y = -1;
  p->slice.z = -1;

}
