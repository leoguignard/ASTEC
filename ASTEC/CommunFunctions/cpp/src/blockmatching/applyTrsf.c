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



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include <chunks.h>
#include <reech4x4.h>
#include <reech-def.h>

#include <bal-stddef.h>
#include <bal-image.h>
#include <bal-transformation.h>
#include <bal-transformation-tools.h>
#include <bal-applyTrsf.h>


/* usages:
   - image_to_be_resampled resampled_image [-trsf %s] [geometry options]
*/


typedef struct {
  char *theim_name;
  char *resim_name;

  char *real_transformation_name;
  char *voxel_transformation_name;

  char *result_real_transformation_name;
  char *result_voxel_transformation_name;

  char *template_image_name; /* template */  
  
  bal_integerPoint dim;
  bal_floatPoint voxel;

  int resize;

  enumTransformationInterpolation interpolation;

  ImageType type;

  int print_time;

} local_parameter;





static char *program = NULL;

static char *usage = "%s %s [-transformation |-trsf %s]\n\
 [-voxel-transformation |-voxel-trsf %s]\n\
 [-template %s] [-dim %d %d [%d]] [-voxel %lf %lf [%lf]]\n\
 [-resize] [-iso %lf]\n\
 [-result-transformation|-res-trsf %s]\n\
 [-result-voxel-transformation|-res-voxel-trsf %s]\n\
 [-nearest|-linear]\n\
 [-parallel|-no-parallel] [-max-chunks %d]\n\
 [-time] [-notime]\n\
 [-o|-b|-bytes %d [-r|-f] [-s]]\n\
 [-verbose|-v | -no-verbose|-nv] [-help]\n";

 

static char *detail = "\
-transformation|-trsf %s # transformation to be applied\n\
  in 'real' coordinates. Goes from 'result image' to 'input image'.\n\
  If indicated, '-voxel-transformation' is ignored.\n\
-voxel-transformation|-voxel-trsf %s #  transformation to be applied\n\
  in 'voxel' coordinates.\n\
-template %s         # template image for the geometry\n\
                       of the output image\n\
-dim %d %d [%d]      # output image dimensions\n\
-voxel %lf %lf [%lf] # voxel sizes of the output image\n\
-resize              # if output image dimensions are given\n\
  output voxel size is computed so that the ouput field of view\n\
  correspond to the input's one (only when no transformation is given)\n\
  if  output image voxel sizes are given, output image dimensions\n\
  are computed so that the ouput field of view correspond to the input's\n\
  one (only when no transformation is given)\n\
-iso %lf # short for '-resize -voxel %lf %lf %lf'\n\
-result-transformation|-res-trsf %s # applied transformation\n\
  in 'real' coordinates. Useful to applied the same transformation to an\n\
  other or to combine transformations.\n\
-result-voxel-transformation|-res-voxel-trsf %s # applied transformation\n\
  in 'voxel' coordinates.\n\
-nearest             # 0 order interpolation, use the value of the nearest\n\
  neighbor\n\
-linear              # 1st order (ie bi- or tri-linear)\n";



static int _verbose_ = 1;
static int _debug_ = 0;


static void _ErrorParse( char *str, int flag );
static void _Parse( int argc, char *argv[], local_parameter *p );
static void _InitParam( local_parameter *p );
static double _GetTime();
static double _GetClock();
static char *_BaseName( char *p );




int main(int argc, char *argv[])
{
  local_parameter p;
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

  /* parsing parameters 
   */
  _InitParam( &p );

  _Parse( argc, argv, &p );

  if ( p.theim_name == NULL ) {
    _ErrorParse( "no input image name", 0 );
    exit( -1 );
  }

  if ( p.resim_name == NULL ) {
    _ErrorParse( "no result image name", 0 );
    exit( -1 );
  }


  

  /***************************************************
   *
   * 
   *
   ***************************************************/

  if (applyTrsf(
		 p.theim_name,
		 p.resim_name,
		 p.real_transformation_name,
		 p.voxel_transformation_name,
		 p.result_real_transformation_name,
		 p.result_voxel_transformation_name,
		 p.template_image_name,
		 p.dim,
		 p.voxel,
		 p.resize,
		 p.interpolation,
		 p.type,
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



static void _Parse( int argc, char *argv[], local_parameter *p )
{
  int i;
  int status;  
  int f=0, o=0, s=0, r=0;
  int maxchunks;

  program = argv[0];

  for ( i=1; i<argc; i++ ) {

    if ( argv[i][0] == '-' ) {

      if ( strcmp ( argv[i], "-template") == 0
	   || (strcmp ( argv[i], "-t") == 0 && argv[i][2] == '\0')
	   ||  strcmp ( argv[i], "-reference") == 0
	   || (strcmp ( argv[i], "-ref") == 0 && argv[i][4] == '\0')
	   || (strcmp ( argv[i], "-dims") == 0 && argv[i][5] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -template", 0 );
	p->template_image_name = argv[i];
      }
      
      else if ( strcmp ( argv[i], "-transformation") == 0
		|| (strcmp ( argv[i], "-trsf") == 0 && argv[i][5] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -transformation", 0 );
	p->real_transformation_name = argv[i];
      }
      
      else if ( strcmp ( argv[i], "-voxel-transformation") == 0
		|| strcmp ( argv[i], "-voxel-trsf") == 0 ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -voxel-transformation", 0 );
	p->voxel_transformation_name = argv[i];
      }
      
      else if ( strcmp ( argv[i], "-result-transformation" ) == 0
		|| strcmp ( argv[i], "-res-trsf" ) == 0 ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -result-transformation", 0 );
	p->result_real_transformation_name = argv[i];
      }
      else if ( strcmp ( argv[i], "-result-voxel-transformation" ) == 0
		|| strcmp ( argv[i], "-res-voxel-trsf" ) == 0 ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -result-voxel-transformation", 0 );
	p->result_voxel_transformation_name = argv[i];
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
      else if ( strcmp (argv[i], "-voxel" ) == 0 && argv[i][6] == '\0' ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -voxel %f", 0 );
	status = sscanf( argv[i], "%f", &(p->voxel.x) );
	if ( status <= 0 ) _ErrorParse( "parsing -voxel %f", 0 );
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -voxel %f %f", 0 );
	status = sscanf( argv[i], "%f", &(p->voxel.y) );
	if ( status <= 0 ) _ErrorParse( "parsing -voxel %f %f", 0 );
	i ++;
	if ( i >= argc) p->voxel.z = 1;
	else {
	  status = sscanf( argv[i], "%f", &(p->voxel.z) );
	  if ( status <= 0 ) {
	    i--;
	    p->voxel.z = 1;
	  }
	}
      }

      else if ( strcmp (argv[i], "-iso" ) == 0 && argv[i][4] == '\0' ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -iso %f", 0 );
	status = sscanf( argv[i], "%f", &(p->voxel.x) );
	if ( status <= 0 ) _ErrorParse( "parsing -voxel %f", 0 );
	p->voxel.y = p->voxel.z = p->voxel.x;
	p->resize = 1;
      }

      else if ( strcmp ( argv[i], "-resize" ) == 0 ) {
	p->resize = 1;
      }

      else if ( strcmp ( argv[i], "-nearest" ) == 0 ) {
	p->interpolation = NEAREST;
      }
      else if ( strcmp ( argv[i], "-linear" ) == 0 ) {
	p->interpolation = LINEAR;
      }


      else if ( strcmp ( argv[i], "-parallel" ) == 0 ) {
	setMaxChunks( 100 );
      }
      
      else if ( strcmp ( argv[i], "-no-parallel" ) == 0 ) {
	setMaxChunks( 1 );
      }

      else if ( strcmp ( argv[i], "-max-chunks" ) == 0 ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "-max-chunks", 0 );
	status = sscanf( argv[i], "%d", &maxchunks );
	if ( status <= 0 ) _ErrorParse( "-max-chunks", 0 );
	if ( maxchunks >= 1 ) setMaxChunks( maxchunks );
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

      else if ( ( strcmp ( argv[i], "-v" ) == 0 && argv[i][2] == '\0' ) 
		|| strcmp ( argv[i], "-verbose" ) == 0 ) {
	if ( _verbose_ <= 0 )  _verbose_ = 1;
	else _verbose_ ++;
	incrementVerboseInReech4x4();
	incrementVerboseInReechDef();
	BAL_IncrementVerboseInBalImage();
      }

      else if ( ( strcmp ( argv[i], "-nv" ) == 0 && argv[i][3] == '\0' ) 
		|| strcmp ( argv[i], "-no-verbose" ) == 0 ) {
	_verbose_ = 0;
	setVerboseInReech4x4( 0 );
	setVerboseInReechDef( 0 );
	BAL_IncrementVerboseInBalImage( 0 );
      }

      else if ( ( strcmp ( argv[i], "-d" ) == 0 && argv[i][2] == '\0' ) 
		|| ( strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0' )) {
	if ( _debug_ <= 0 )  _debug_ = 1;
	else _debug_ ++;
      }
      
      /* image encoding type
       */

      else if ( strcmp ( argv[i], "-f" ) == 0 && argv[i][2] == '\0' ) {
	f = 1;
      }
      else if ( strcmp ( argv[i], "-r" ) == 0 && argv[i][2] == '\0' ) {
	r = 1;
      }
      else if ( strcmp ( argv[i], "-s" ) == 0 && argv[i][2] == '\0' ) {
	s = 1;
      }
      else if ( strcmp ( argv[i], "-bytes" ) == 0
		|| (strcmp ( argv[i], "-b" ) == 0 && argv[i][2] == '\0')
		|| (strcmp ( argv[i], "-o" ) == 0 && argv[i][2] == '\0') ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -bytes...\n", 0 );
	status = sscanf( argv[i],"%d",&o );
	if ( status <= 0 ) _ErrorParse( "parsing -bytes...\n", 0 );
      }

      else {
	fprintf(stderr,"unknown option: '%s'\n",argv[i]);
      }

    }
    /* ( argv[i][0] == '-' )
     */
    else {
      if ( p->theim_name == NULL ) {
	p->theim_name = argv[i];
      }
      else if ( p->resim_name == NULL ) {
	p->resim_name = argv[i];
      }
      else {
	fprintf(stderr,"too many file names: '%s'\n",argv[i]);
      }
    }

  }
  
   /*--- type de l'image resultat ---*/
  if ( r == 1 ) {
    switch( o ) {
    default :
      fprintf( stderr, "such byte size not handled for floating types\n" );
      break;
    case 0 :
    case 4 :
      p->type = FLOAT; break;
    case 8 :
      p->type = FLOAT; break;
    }
  }
  else {
    switch( o ) {
    default :
      _ErrorParse( "such byte size not handled for integer types\n", 0 );
      break;
    case 0 :
      /* ce 'break' est necessaire, sinon le type par defaut est 'UCHAR'
       */
      break;
    case 1 :
      p->type = ( s == 1 ) ? SCHAR : UCHAR;
      break;
    case 2 :
      p->type = ( s == 1 ) ? SSHORT : USHORT;
      break;
    case 4 :
      p->type = ( s == 1 ) ? SINT : UINT;
      break;
    case 8 :
      if ( s == 1 )
	_ErrorParse( "signed long int not handled yet\n", 0 );
      else 
	p->type = ULINT;
      break;
    }
  }


}



static void _InitParam( local_parameter *p )
{
  p->theim_name = NULL;
  p->resim_name = NULL;

  p->real_transformation_name = NULL;
  p->voxel_transformation_name = NULL;

  p->result_real_transformation_name = NULL;
  p->result_voxel_transformation_name = NULL;

  p->template_image_name = NULL;
  
  p->dim.x = -1;
  p->dim.y = -1;
  p->dim.z = -1;

  p->voxel.x = -1.0;
  p->voxel.y = -1.0;
  p->voxel.z = -1.0;

  p->resize = 0;

  p->interpolation = LINEAR;

  p->type = TYPE_UNKNOWN;

  p->print_time = 0;
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
