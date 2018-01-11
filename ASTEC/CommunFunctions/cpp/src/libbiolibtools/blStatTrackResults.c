/*************************************************************************
 * blStatTrackResults.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2014, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Lun  3 mar 2014 16:27:28 CET
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



#include <histogram.h>

#include <bal-image.h>

#include <bal-biolib-tools.h>






typedef struct local_parameter {

  /* file names
   */
  char *thetrack_name;
  char *resname;

  char *template_image_name;

  /* voxel size
   */
  bal_doublePoint voxel;
  double deltat;

  enumHistogramFile xxxlab;

} local_parameter;





/*------- Definition des fonctions statiques ----------*/
static void _ErrorParse( char *str, int flag );
static void _Parse( int argc, char *argv[], local_parameter *p );
static void _InitParam( local_parameter *par );




static char *program = NULL;

static char *usage = "%s -scilab %s \n\
 [-template %s]\n\
 [-voxel | -pixel | -vs %f %f [%f] | [-vx %f] [-vy %f] [-vz %f] ]\n\
 [-deltat | -dt %lf]\n\
 [-v] [-help]";

static char *detail = "\
-template %s         # template image for the dimensions\n\
                       of the output image\n\
-voxel %lf %lf [%lf]    # image voxel sizes\n\
-deltat %lf          # time step\n\
 -v : mode verbose\n\
\n";











int main( int argc, char *argv[] )
{
  local_parameter p;
  bal_blTrackList list;
  bal_doublePoint voxelsize;
  bal_image template;

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
  

  /* geometry stuff
   */
  voxelsize.x = voxelsize.y = voxelsize.z = 0.0;
  
  if ( p.template_image_name != NULL ) {
    if ( BAL_ReadImage( &template, p.template_image_name, 0 ) != 1 ) {
      _ErrorParse( "can not read template image\n", 0 );
    }
    voxelsize.x = template.vx;
    voxelsize.y = template.vy;
    voxelsize.z = template.vz;
    BAL_FreeImage( &template );
  }
  else if ( p.voxel.x > 0 &&  p.voxel.y > 0 ) {
    voxelsize.x = p.voxel.x;
    voxelsize.y = p.voxel.y;
    voxelsize.z = p.voxel.z;
    if ( voxelsize.z <= 0.0 ) voxelsize.z = 1.0;
  }
  else {
    voxelsize.x = voxelsize.y = voxelsize.z = 1.0;
  }



  
  BAL_InitBlTrackList( &list );
  if ( p.thetrack_name == (char*)NULL ) {
    _ErrorParse( "no detection file\n", 0 );
  }

  if ( BAL_ReadBlTrackList( &list, p.thetrack_name ) != 0 ) {
    _ErrorParse( "unable to read track file\n", 0 );
  }

  BAL_ComputeTrackListProperties( &list, voxelsize.x, voxelsize.y, voxelsize.z, p.deltat );

  switch( p.xxxlab ) {
  default :
    BAL_FreeBlTrackList( &list );
    _ErrorParse( "such output type not handled yet\n", 0 );
    break;  case _SCILAB_ :
    if ( BAL_PrintAllStatsScilab( p.resname, &list ) != 1 ) {
      BAL_FreeBlTrackList( &list );
      _ErrorParse( "unable to write track file\n", 0 );
    }
    break;
  }


  BAL_FreeBlTrackList( &list );
  return( 1 );
}













static void _Parse( int argc, char *argv[], local_parameter *p )
{
  int i, status;
  int inputisread = 0;
  int outputisread = 0;

  program = argv[0];
	
  for ( i=1; i<argc; i++ ) {
  
    if ( argv[i][0] == '-' ) {

      if ( strcmp ( argv[i], "-template") == 0
	   || (strcmp ( argv[i], "-t") == 0 && argv[i][2] == '\0')
	   || (strcmp ( argv[i], "-dims") == 0 && argv[i][5] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -template\n", 0 );
	p->template_image_name = argv[i];
      }

      else if ( strcmp ( argv[i], "-voxel" ) == 0
		|| strcmp ( argv[i], "-pixel" ) == 0
		|| (strcmp ( argv[i], "-vs" ) == 0  && argv[i][3] == '\0') ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -voxel..\n", 0 );
	status = sscanf( argv[i],"%lf",&(p->voxel.x) );
	if ( status <= 0 ) _ErrorParse( "parsing -voxel...\n", 0 );
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -voxel..\n", 0 );
	status = sscanf( argv[i],"%lf",&(p->voxel.y) );
	if ( status <= 0 ) _ErrorParse( "parsing -voxel...\n", 0 );
	i += 1;
	if ( i < argc ) {
	  status = sscanf( argv[i],"%lf",&(p->voxel.z) );
	  if ( status <= 0 ) {
	    i--;
	  }
	}
      }

      else if ( strcmp ( argv[i], "-vx" ) == 0   && argv[i][3] == '\0' ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -vx...\n", 0 );
	status = sscanf( argv[i], "%lf", &(p->voxel.x) );
	if ( status <= 0 ) _ErrorParse( "parsing -vx...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-vy" ) == 0   && argv[i][3] == '\0' ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -vy...\n", 0 );
	status = sscanf( argv[i], "%lf", &(p->voxel.y) );
	if ( status <= 0 ) _ErrorParse( "parsing -vy...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-vz" ) == 0   && argv[i][3] == '\0' ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -vz...\n", 0 );
	status = sscanf( argv[i], "%lf", &(p->voxel.z) );
	if ( status <= 0 ) _ErrorParse( "parsing -vz...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-deltat" ) == 0 
		|| (strcmp ( argv[i], "-dt" ) == 0   && argv[i][3] == '\0') ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -deltat...\n", 0 );
	status = sscanf( argv[i], "%lf", &(p->deltat) );
	if ( status <= 0 ) _ErrorParse( "parsing -deltat...\n", 0 );
      }



      /* output type
       */
      
      else if ( strcmp ( argv[i], "-matlab" ) == 0 ) {
	p->xxxlab = _MATLAB_;
      }

      else if ( strcmp ( argv[i], "-scilab" ) == 0 ) {
	p->xxxlab = _SCILAB_;
      }


      /* general arguments
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
	p->resname = argv[i];
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
  p->resname = (char*)NULL;

  p->template_image_name = (char*)NULL;
  
  p->voxel.x = -1.0;
  p->voxel.y = -1.0;
  p->voxel.z = -1.0;
  p->deltat = 1;
}
