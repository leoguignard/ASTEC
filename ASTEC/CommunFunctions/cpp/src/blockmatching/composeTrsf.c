/*************************************************************************
 * composeTrsf.c -
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
#include <sys/time.h>

#include <bal-stddef.h>
#include <bal-transformation.h>
#include <bal-transformation-tools.h>

#include <bal-composeTrsf.h>

typedef struct {
  char *restrsf_name;

  char *template_image_name; /* template */  
  
  bal_integerPoint dim;
  bal_doublePoint voxel;

  /*
  int first_trsf;
  int last_trsf;
  */

  int print_time;

} local_parameter;








static char *program = NULL;

static char *usage = "[-res output-transformation] [-trsfs %s ... %s]\n\
 [-template %s] [-dim %d %d [%d]] [-voxel %lf %lf [%lf]]\n\
 [-time] [-notime]\n";

static char *detail = "\
-res|-result %s      # output-transformation\n\
-template %s         # template image for geometry\n\
-dim %d %d [%d]      # template image dimensions\n\
-voxel %lf %lf [%lf] # voxel sizes of the template image\n\
\n\
Transformations are composed in the order they are given.\n\
The line '-trsfs T1 T2 ... TN' assumes that transformation\n\
Ti goes from image I(i+1) to I(i) [then allows to resample\n\
I(i) in the same frame than I(i+1)], ie  Ti = T_{I(i) <- I(i+1)}.\n\
The resulting transformation will goes from I(N+1) to I(1)\n\
[then allows to resample I(1) in the same frame than I(N+1)]\n";



static int _verbose_ = 1;
static int _debug_ = 0;


static void _ErrorParse( char *str, int flag );
static void _Parse( int argc, char *argv[], local_parameter *p, int *is_a_trsf );
static void _PrintArgs( FILE *f, int argc, char *argv[], int *is_a_trsf );
static void _InitParam( local_parameter *p );
static double _GetTime();
static double _GetClock();
static char *_BaseName( char *p );





int main(int argc, char *argv[])
{
  local_parameter p;
  int i;
  int *is_a_trsf = (int*)NULL;
  
  double time_init = _GetTime();
  double time_exit;
  double clock_init = _GetClock();
  double clock_exit;


  /***************************************************
   *
   * parsing parameters
   *
   ***************************************************/
  program = argv[0];
  
  
  /* no arguments
   */
  if ( argc == 1 ) _ErrorParse( NULL, 0 );

  
  /* displaying help is required
   */
  i = 1;
  while ( i < argc ) {
    if ( ( strcmp ( argv[i], "-help") == 0 ) 
	 || ( strcmp ( argv[i], "-h") == 0 && argv[i][2] == '\0' ) 
	 || ( strcmp ( argv[i], "--help") == 0 ) 
	 || ( strcmp ( argv[i], "--h") == 0 && argv[i][3] == '\0' ) ) {
      _ErrorParse( NULL, 1 );
    }
    i++;
  }

  is_a_trsf = (int*)malloc( argc*sizeof(int) );
  if ( is_a_trsf == (int*) NULL ) {
    _ErrorParse( "error when allocating auxiliary array", 0 );
  }
  for ( i=0; i<argc; i++ ) is_a_trsf[i] = 0;
  
  _InitParam( &p );
  _Parse( argc, argv, &p, is_a_trsf );
  
  
  /* check the transformations to be composed
   */

  if ( _debug_ ) _PrintArgs( stderr, argc, argv, is_a_trsf );

  if(composeTrsf(
        p.restrsf_name,
        p.template_image_name,
        p.dim,
        p.voxel,
	/*
	  p.first_trsf,
	  p.last_trsf,
	*/
	argv,
	argc,
	is_a_trsf,
	_debug_,
	_verbose_
    ))
  {
    fprintf( stderr, "%s: Failure.\n",program);
  }
  
  time_exit = _GetTime();
  clock_exit = _GetClock();

  if ( p.print_time ) { 
    fprintf( stderr, "%s: elapsed (real) time = %f\n", _BaseName(program), time_exit - time_init );
    fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_init );
    fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_init)/(time_exit - time_init) );
  }
  
  if ( is_a_trsf != (int*)NULL )
    free( is_a_trsf );

  exit( 0 );
}




/***************************************************
 *
 * 
 *
 ***************************************************/





static void _ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage: %s %s\n",_BaseName(program), usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s\n",detail);
  if ( str != NULL ) (void)fprintf(stderr,"Error: %s\n",str);
  exit( -1 );
}

static void _Parse( int argc, char *argv[], local_parameter *p, int *is_a_trsf )
{
  int i;
  int status;  

  program = argv[0];

  for ( i=1; i<argc; i++ ) {

    if ( argv[i][0] == '-' ) {

      if ( strcmp ( argv[i], "-result") == 0
	   || (strcmp ( argv[i], "-res") == 0 && argv[i][4] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -result", 0 );
	p->restrsf_name = argv[i];
      }
      
      else if ( strcmp ( argv[i], "-template") == 0
		|| (strcmp ( argv[i], "-t") == 0 && argv[i][2] == '\0')
		|| (strcmp ( argv[i], "-dims") == 0 && argv[i][5] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -reference", 0 );
	p->template_image_name = argv[i];
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
      else if ( strcmp (argv[i], "-voxel" ) == 0 ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -voxel %lf", 0 );
	status = sscanf( argv[i], "%lf", &(p->voxel.x) );
	if ( status <= 0 ) _ErrorParse( "parsing -voxel %lf", 0 );
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -voxel %lf %lf", 0 );
	status = sscanf( argv[i], "%lf", &(p->voxel.y) );
	if ( status <= 0 ) _ErrorParse( "parsing -voxel %lf %lf", 0 );
	i ++;
	if ( i >= argc) p->voxel.z = 1;
	else {
	  status = sscanf( argv[i], "%lf", &(p->voxel.z) );
	  if ( status <= 0 ) {
	    i--;
	    p->voxel.z = 1;
	  }
	}
      }

      else if ( strcmp ( argv[i], "-transformations" ) == 0 
		|| strcmp ( argv[i], "-trsfs" ) == 0 ) {
	i++;
	while ( i<argc && argv[i][0] != '-' ) {
	  is_a_trsf[i] = 1;
	  i++;
	}
	i--;
      }
      
      else if ( (strcmp ( argv[i], "-time" ) == 0 && argv[i][5] == '\0') ) {
	p->print_time = 1;
      }
      else if ( (strcmp ( argv[i], "-notime" ) == 0 && argv[i][7] == '\0') 
		|| (strcmp ( argv[i], "-no-time" ) == 0 && argv[i][8] == '\0') ) {
	p->print_time = 0;
      }

      else if ( strcmp ( argv[i], "--help" ) == 0 
		|| ( strcmp ( argv[i], "-help" ) == 0 && argv[i][5] == '\0' )
		|| ( strcmp ( argv[i], "--h" ) == 0 && argv[i][3] == '\0' )
		|| ( strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0' ) ) {
	_ErrorParse( NULL, 1 );
      }

      else if ( ( strcmp ( argv[i], "-d" ) == 0 && argv[i][2] == '\0' ) 
		|| ( strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0' )) {
	if ( _debug_ <= 0 )  _debug_ = 1;
	else _debug_ ++;
      }

      else {
	fprintf(stderr,"unknown option: '%s'\n",argv[i]);
      }

    }
    /* ( argv[i][0] == '-' )
     */
    else {
      fprintf(stderr,"too many file names: '%s'\n",argv[i]);
    }

  }
  
}



static void _PrintArgs( FILE *f, int argc, char *argv[], int *is_a_trsf )
{
  int i;
  fprintf( f, "\n" );
  fprintf( f, " flag # 1 if it a transformation to be composed\n" );
  fprintf( f, " ----------------------------------------------\n" );
  for ( i=0; i<argc; i++ )
    fprintf( f, " %d # '%s'\n", is_a_trsf[i], argv[i] );
  fprintf( f, "\n" );
}



static void _InitParam( local_parameter *p )
{
  p->restrsf_name = NULL;

  p->template_image_name = NULL;
  
  p->dim.x = -1;
  p->dim.y = -1;
  p->dim.z = -1;

  p->voxel.x = -1.0;
  p->voxel.y = -1.0;
  p->voxel.z = -1.0;

  /*
  p->first_trsf = -1;
  p->last_trsf = -1;
  */

  p->print_time = 1;
}





static double _GetTime() 
{
  struct timeval tv;
  gettimeofday(&tv, (void *)0);
  return ( (double) tv.tv_sec + tv.tv_usec*1e-6 );
}

static double _GetClock() 
{
  return ( (double) clock() / (double)CLOCKS_PER_SEC );
}





static char *_BaseName( char *p )
{
  int l;
  if ( p == (char*)NULL ) return( (char*)NULL );
  l = strlen( p ) - 1;
  while ( l >= 0 && p[l] != '/' ) l--;
  if ( l < 0 ) l = 0;
  if ( p[l] == '/' ) l++;
  return( &(p[l]) );
}

