/*************************************************************************
 * build-pyramid-image.c -
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

#include <sys/time.h>
#include <time.h>


#include <chunks.h>

#include <bal-image.h>
#include <bal-image-tools.h>
#include <bal-pyramid.h>

/*
static int _debug_ = 0;
*/

/*------- Definition des fonctions statiques ----------*/

static char *program = NULL;

/*---------------- longueur maximale de ligne --------------------------------*/
/*----------------------------------------------------------------------------*/
static char *usage = "-image-prefix %s -image-suffix %s\n\
 -transformation-prefix|-trsf-prefix %s\n\
 [-normalisation|-norma|-rescale|-no-normalisation|-no-norma|-no-rescale]\n\
 [-pyramid-lowest-level | -py-ll %d] [-pyramid-highest-level | -py-hl %d]\n\
 [-pyramid-gaussian-filtering | -py-gf]\n\
 [-gaussian-filter-type|-filter-type deriche|fidrich|young-1995|young-2002|...\n\
 ...|gabor-young-2002|convolution]\n\
 [-parallel|-no-parallel] [-max-chunks %d]\n\
 [-parallel-scheduling|-ps default|static|dynamic-one|dynamic|guided]\n\
 [-verbose|-v] [-no-verbose|-nv] [-time] [-notime]\n\
 [-h|-help|--h|--help]\n";

static char *detail = "\
 -image-prefix %s # \n\
 -image-suffix %s # \n\
 -transformation-prefix|-trsf-prefix %s # \n\
[-normalisation|-norma|-rescale] # input images are normalized on one byte\n\
  before matching (this may be the default behavior)\n\
[-no-normalisation|-no-norma|-no-rescale] # input images are not normalized on\n\
   one byte before matching\n\
### pyramid building ###\n\
[-pyramid-lowest-level | -py-ll %d]    # pyramid lowest level (0 = original dimension)\n\
[-pyramid-highest-level | -py-hl %d]   # pyramid highest level\n\
  default is 3: it corresponds to 32x32x32 for an original 256x256x256 image\n\
[-pyramid-gaussian-filtering | -py-gf] # before subsampling, the images are filtered\n\
### filter type ###\n\
[-gaussian-filter-type|-filter-type deriche|fidrich|young-1995|young-2002|...\n\
  ...|gabor-young-2002|convolution] # type of filter for image/vector field smoothing\n\
[-verbose] [-no-verbose]\n\
[-h|-help|--h|--help]\n\
\n";


#define STRLENGTH 1024 


typedef struct local_parameter {
	
  /* file names
     - images
     - transformations
  */
  char *input_image;
  char output_prefix_image[STRLENGTH];
  char output_suffix_image[STRLENGTH];
  char output_prefix_trsf[STRLENGTH];
  
  int normalisation;

  int pyramid_lowest_level;
  int pyramid_highest_level;
  int pyramid_gaussian_filtering;

  int verbose;
  int print_time;

} local_parameter;



#define FILENAMELENGTH 128



static void _ErrorParse( char *str, int flag );
static void _Parse( int argc, char *argv[], local_parameter *p );
static void _InitParam( local_parameter *par );
static double _GetTime();
static double _GetClock();
static char *_BaseName( char *p );








int main(int argc, char *argv[])
{
  local_parameter p;

  bal_image theImage;

  double time_init = _GetTime();
  double time_exit;
  double clock_init = _GetClock();
  double clock_exit;



  /* parsing parameters 
   */
  _InitParam( &p );
  _Parse( argc, argv, &p );



  if ( p.input_image == (char*)NULL ) {
    _ErrorParse( "no input image", 0 );
  }

  BAL_InitImage( &theImage, NULL, 0, 0, 0, 0, UCHAR );
  if ( BAL_ReadImage( &theImage, p.input_image, p.normalisation ) != 1 ) {
    fprintf( stderr, "%s: can not read '%s'\n", _BaseName(argv[0]), p.input_image );
    exit( -1 );
  }


  if ( BAL_BuildPyramidImage( &theImage,
			      p.output_prefix_image,
			      p.output_suffix_image,
			      p.output_prefix_trsf,
			      p.pyramid_lowest_level,
			      p.pyramid_highest_level,
			      p.pyramid_gaussian_filtering ) != 1 ) {       
    BAL_FreeImage( &theImage );
    fprintf( stderr, "%s: unable to build pyramid\n", _BaseName(argv[0]) );
    exit( -1 );
  }

  BAL_FreeImage( &theImage );


 
  time_exit = _GetTime();
  clock_exit = _GetClock();

  if ( p.print_time ) { 
    fprintf( stderr, "%s: elapsed (real) time = %f\n", _BaseName(program), time_exit - time_init );
    fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_init );
    fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_init)/(time_exit - time_init) );
  }

  return( 0 );
}







/************************************************************
 *
 *
 *
 ************************************************************/


static void _ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage: %s %s\n",_BaseName(program), usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s\n",detail);
  if ( str != NULL ) (void)fprintf(stderr,"Error: %s\n",str);
  exit( -1 );
}



/************************************************************
 *
 *
 * 
 ************************************************************/

static void _Parse( int argc, char *argv[], local_parameter *p )
{
  int i;
  int status;
  int maxchunks;

  program = argv[0];
  


  /* reading the transformation type
   */

  for ( i=1; i<argc; i++ ) {

    if ( argv[i][0] == '-' ) {

      /* names
       */
      if ( strcmp ( argv[i], "-image-prefix" ) == 0 ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "-image-prefix", 0 );
	(void)strncpy( p->output_prefix_image, argv[i], STRLENGTH );
      }
      else if ( strcmp ( argv[i], "-image-suffix" ) == 0 ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "-image-suffix", 0 );
	(void)strncpy( p->output_suffix_image, argv[i], STRLENGTH );
      }
      else if ( strcmp ( argv[i], "-transformation-prefix" ) == 0 
	   || strcmp ( argv[i], "-trsf-prefix" ) == 0 ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "-transformation-prefix", 0 );
	(void)strncpy( p->output_prefix_trsf, argv[i], STRLENGTH );
      }
      
      /* pre-processing before matching
	 input images normalisation
      */
      else if ( strcmp ( argv[i], "-normalisation") == 0
		|| (strcmp ( argv[i], "-norma") == 0 && argv[i][6] == '\0') 
		|| (strcmp ( argv[i], "-rescale") == 0 && argv[i][8] == '\0') ) {
	p->normalisation = 1;
      }
      else if ( strcmp ( argv[i], "-no-normalisation") == 0 
		|| (strcmp ( argv[i], "-no-norma") == 0 && argv[i][9] == '\0') 
		|| (strcmp ( argv[i], "-no-rescale") == 0 && argv[i][11] == '\0') ) {
	p->normalisation = 0;
      }
      
      
      /* pyramid building
       */
      else if ( strcmp ( argv[i], "-pyramid-lowest-level" ) == 0
		|| (strcmp( argv[i], "-py-ll") == 0 && argv[i][6] == '\0') ) { 
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -pyramid-lowest-level", 0 );
	status = sscanf( argv[i], "%d", &(p->pyramid_lowest_level) );
	if ( status <= 0 ) _ErrorParse( "parsing -pyramid-lowest-level", 0 );
      }
      else if ( strcmp ( argv[i], "-pyramid-highest-level" ) == 0
		|| (strcmp( argv[i], "-py-hl") == 0 && argv[i][6] == '\0') ) { 
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -pyramid-highest-level", 0 );
	status = sscanf( argv[i], "%d", &(p->pyramid_highest_level) );
	if ( status <= 0 ) _ErrorParse( "parsing -pyramid-highest-level", 0 );
      }
      else if ( strcmp ( argv[i], "-pyramid-gaussian-filtering" ) == 0
		|| (strcmp( argv[i], "-py-gf") == 0 && argv[i][6] == '\0') ) { 
	p->pyramid_gaussian_filtering = 1;
      }
      
      
      /* filter type for image smoothing
       */
      else if ( strcmp ( argv[i], "-gaussian-filter-type" ) == 0 
		|| strcmp ( argv[i], "-filter-type" ) == 0 ) {
	i++;
	if ( strcmp ( argv[i], "deriche" ) == 0 ) {
	  BAL_SetFilterType( GAUSSIAN_DERICHE );
	}
	else if ( strcmp ( argv[i], "fidrich" ) == 0 ) {
	  BAL_SetFilterType( GAUSSIAN_FIDRICH );
	}
	else if ( strcmp ( argv[i], "young-1995" ) == 0 ) {
	  BAL_SetFilterType( GAUSSIAN_YOUNG_1995 );
	}
	else if ( strcmp ( argv[i], "young-2002" ) == 0 ) {
	  BAL_SetFilterType( GAUSSIAN_YOUNG_2002 );
	}
	else if ( strcmp ( argv[i], "gabor-young-2002" ) == 0 ) {
	  BAL_SetFilterType( GABOR_YOUNG_2002 );
	}
	else if ( strcmp ( argv[i], "convolution" ) == 0  ) {
	  BAL_SetFilterType( GAUSSIAN_CONVOLUTION );
	}
	else {
	  fprintf( stderr, "unknown filter type: '%s'\n", argv[i] );
	  _ErrorParse( "-gaussian-filter-type", 0 );
	}
      }
      
      
      /* parallelism
       */
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
      
      else if ( strcmp ( argv[i], "-parallel-scheduling" ) == 0 || 
		( strcmp ( argv[i], "-ps" ) == 0 && argv[i][3] == '\0') ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "-parallel-scheduling", 0 );
	if ( strcmp ( argv[i], "default" ) == 0 ) {
	  setOpenMPScheduling( _DEFAULT_SCHEDULING_ );
	}
	else if ( strcmp ( argv[i], "static" ) == 0 ) {
	  setOpenMPScheduling( _STATIC_SCHEDULING_ );
	}
	else if ( strcmp ( argv[i], "dynamic-one" ) == 0 ) {
	  setOpenMPScheduling( _DYNAMIC_ONE_SCHEDULING_ );
	}
	else if ( strcmp ( argv[i], "dynamic" ) == 0 ) {
	  setOpenMPScheduling( _DYNAMIC_SCHEDULING_ );
	}
	else if ( strcmp ( argv[i], "guided" ) == 0 ) {
	  setOpenMPScheduling( _GUIDED_SCHEDULING_ );
	}
	else {
	  fprintf( stderr, "unknown scheduling type: '%s'\n", argv[i] );
	  _ErrorParse( "-parallel-scheduling", 0 );
	}
      }
      
      
      
      /* misc writing stuff
       */
      
      else if ( (strcmp ( argv[i], "-time" ) == 0 && argv[i][5] == '\0') ) {
	p->print_time = 1;
      }
      else if ( (strcmp ( argv[i], "-notime" ) == 0 && argv[i][7] == '\0') 
	      || (strcmp ( argv[i], "-no-time" ) == 0 && argv[i][8] == '\0') ) {
	p->print_time = 0;
      }
      
      
      
      else if ( (strcmp ( argv[i], "-verbose") == 0 && argv[i][8] == '\0')
		|| (strcmp ( argv[i], "-v") == 0 && argv[i][2] == '\0') ) {
	if ( p->verbose <= 0 ) p->verbose = 1;
	else p->verbose ++;
	
	BAL_IncrementVerboseInBalImageTools(  );
	BAL_IncrementVerboseInBalImage(  );
	
      }
      else if ( (strcmp ( argv[i], "-no-verbose") == 0) 
		|| (strcmp ( argv[i], "-nv") == 0 && argv[i][3] == '\0') ) {
	p->verbose = 0;
	
	BAL_DecrementVerboseInBalImageTools(  );
	BAL_DecrementVerboseInBalImage(  );
	
      }
      
      else if ( ( strcmp ( argv[i], "--help" ) == 0 && argv[i][6] == '\0' )
		|| ( strcmp ( argv[i], "-help" ) == 0 && argv[i][5] == '\0' ) ) {
	_ErrorParse( NULL, 1 );
      }
      
      else if ( ( strcmp ( argv[i], "--h" ) == 0 && argv[i][3] == '\0' )
		|| ( strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0' ) ) {
	_ErrorParse( NULL, 0 );
      }
      
      else {
	fprintf(stderr,"unknown option: '%s'\n",argv[i]);
      }

    }	
    /* ( argv[i][0] == '-' )
     */
    else {
      if ( p->input_image == NULL ) {
	p->input_image = argv[i];
      }
      else {
	fprintf(stderr,"too many file names: '%s'\n",argv[i]);
      }
    }

  }
  

}





static void _InitParam( local_parameter *p ) 
{
  /* file names
     - images
     - transformations
  */

  p->input_image = (char*)NULL;
  p->output_prefix_image[0] = '\0';
  p->output_suffix_image[0] = '\0';
  p->output_prefix_trsf[0] = '\0';

  p->normalisation = 0;

  p->pyramid_lowest_level = 0;
  p->pyramid_highest_level = -1;
  p->pyramid_gaussian_filtering = 0;

  p->verbose = 1;
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


