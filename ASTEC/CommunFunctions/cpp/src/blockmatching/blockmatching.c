/*************************************************************************
 * blockmatching.c -
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
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include <chunks.h>

#include <bal-behavior.h>
#include <bal-blockmatching-param.h>
#include <bal-blockmatching-param-tools.h>
#include <bal-blockmatching.h>
#include <bal-transformation-tools.h>
#include <bal-field-tools.h>

#include <bal-blockmatching_fct.h>

/* for verbose 
 */

#include <recline.h>
#include <reech-def.h>
#include <reech4x4.h>

#include <bal-block-tools.h>
#include <bal-image-tools.h>
#include <bal-lineartrsf.h>
#include <bal-pyramid.h>
#include <bal-vectorfield.h>




/*
  #define max(a,b) ((a)>(b) ? (a) : (b))
  #define min(a,b) ((a)<(b) ? (a) : (b))
*/


static int _debug_ = 0;

/*------- Definition des fonctions statiques ----------*/

static char *program = NULL;

/*---------------- longueur maximale de ligne --------------------------------*/
/*----------------------------------------------------------------------------*/
static char *usage = "-reference|-ref %s -floating|-flo %s -result|-res %s\n\
 [-initial-transformation|-init-trsf %s]\n\
 [-initial-voxel-transformation|-init-voxel-trsf %s]\n\
 [-initial-result-transformation|-init-res-trsf %s]\n\
 [-initial-result-voxel-transformation|-init-res-voxel-trsf %s]\n\
 [-result-transformation|-res-trsf %s]\n\
 [-result-voxel-transformation|-res-voxel-trsf %s]\n";

static char *options = "\
 [-normalisation|-norma|-rescale|-no-normalisation|-no-norma|-no-rescale]\n\
 [-pyramid-lowest-level | -py-ll %d] [-pyramid-highest-level | -py-hl %d]\n\
 [-pyramid-gaussian-filtering | -py-gf]\n\
 [-block-size|-bl-size %d %d %d] [-block-spacing|-bl-space %d %d %d]\n\
 [-block-border|-bl-border %d %d %d]\n\
 [-floating-low-threshold | -flo-lt %d]\n\
 [-floating-high-threshold | -flo-ht %d]\n\
 [-floating-removed-fraction | -flo-rf %f]\n\
 [-reference-low-threshold | -ref-lt %d]\n\
 [-reference-high-threshold | -ref-ht %d]\n\
 [-reference-removed-fraction | -ref-rf %f]\n\
 [-floating-selection-fraction[-ll|-hl] | -flo-frac[-ll|-hl] %lf]\n\
 [-search-neighborhood-half-size | -se-hsize %d %d %d]\n\
 [-search-neighborhood-step | -se-step %d %d %d]\n\
 [-similarity-measure | -similarity | -si [cc]]\n\
 [-similarity-measure-threshold | -si-th %lf]\n\
 [-transformation-type|-transformation|-trsf-type %s]\n\
 [-elastic-regularization-sigma[-ll|-hl] | -elastic-sigma[-ll|-hl]  %lf %lf %lf]\n\
 [-estimator-type|-estimator|-es-type wlts|lts|wls|ls]\n\
 [-lts-cut|-lts-fraction %lf] [-lts-deviation %f] [-lts-iterations %d]\n\
 [-fluid-sigma|-lts-sigma[-ll|-hl] %lf %lf %lf]\n\
 [-max-iteration[-ll|-hl]|-max-iter[-ll|-hl] %d] [-corner-ending-condition|-rms]\n\
 [-gaussian-filter-type|-filter-type deriche|fidrich|young-1995|young-2002|...\n\
 ...|gabor-young-2002|convolution]\n\
 [-parallel|-no-parallel] [-max-chunks %d]\n\
 [-parallel-scheduling|-ps default|static|dynamic-one|dynamic|guided]\n\
 [-print-parameters]\n\
 [-default-filenames|-df] [-no-default-filenames|-ndf]\n\
 [-command-line %s] [-logfile %s]\n\
 [-vischeck] [-write_def]\n\
 [-verbose|-v] [-no-verbose|-nv] [-time] [-notime]\n\
 [-debug|-no-debug] [-trace|-no-trace]\n\
 [-h|-help|--h|--help]\n";

static char *detail = "\
### File names ###\n\
-reference|-ref %s  # name of the reference image (still image)\n\
-floating|-flo %s   # name of the image to be registered (floating image)\n\
  blocks are defined on this image\n\
-result|res %s      # name of the result image (default is 'res.inr.gz')\n\
  this is the floating image resampled with the output transformation\n\
  in the same geometry than the reference image\n\
[-initial-transformation|-init-trsf %s] # name of the initial transformation\n\
  in 'real' coordinates. Goes from 'reference' to 'floating', ie allows to \n\
  resample 'floating' in the geometry of 'reference' before registration.\n\
  Comes then to register 'floating o init' with 'reference'.\n\
  If indicated, '-initial-voxel-transformation' is ignored.\n\
[-initial-voxel-transformation|-init-voxel-trsf %s] # name of the initial\n\
  transformation in 'voxel' coordinates.\n\
[-initial-result-transformation|-init-res-trsf %s] # name of the\n\
  initialization of the result transformation in 'real' coordinates.\n\
  May be used when the registration has only be conducted at high scales\n\
  and one want to continue it at lower scales.\n\
  Comes then to register 'floating o init' with 'reference'.\n\
  If indicated, '-initial-result-voxel-transformation' is ignored.\n\
[-initial-result-voxel-transformation|-init-res-voxel-trsf %s] # name of the\n\
  initialization of the result transformation in 'voxel' coordinates.\n\
[-result-transformation|-res-trsf %s] # name of the result transformation\n\
  in 'real' coordinates. Goes from 'reference' to 'floating', ie allows to \n\
  resample 'floating' in the geometry of 'reference'.\n\
  If indicated, '-result-voxel-transformation' is ignored.\n\
[-result-voxel-transformation|-res-voxel-trsf %s] # name of the result\n\
  transformation in 'voxel' coordinates.\n\
### pre-processing ###\n\
[-normalisation|-norma|-rescale] # input images are normalized on one byte\n\
  before matching (this may be the default behavior)\n\
[-no-normalisation|-no-norma|-no-rescale] # input images are not normalized on\n\
  one byte before matching\n\
### pyramid building ###\n\
[-pyramid-lowest-level | -py-ll %d]    # pyramid lowest level\n\
  (0 = original dimension)\n\
[-pyramid-highest-level | -py-hl %d]   # pyramid highest level\n\
  default is 3: it corresponds to 32x32x32 for an original 256x256x256 image\n\
[-pyramid-gaussian-filtering | -py-gf] # before subsampling, the images \n\
  are filtered (ie smoothed) by a gaussian kernel.\n\
### block geometry (floating image) ###\n\
-block-size|-bl-size %d %d %d       # size of the block along X, Y, Z\n\
-block-spacing|-bl-space %d %d %d   # block spacing in the floating image\n\
-block-border|-bl-border %d %d %d   # block borders: to be added twice at\n\
  each dimension for statistics computation\n\
### block selection ###\n\
[-floating-low-threshold | -flo-lt %d]     # values <= low threshold are not\n\
  considered\n\
[-floating-high-threshold | -flo-ht %d]    # values >= high threshold are not\n\
  considered\n\
[-floating-removed-fraction | -flo-rf %f]  # maximal fraction of removed points\n\
  because of the threshold. If too many points are removed, the block is\n\
  discarded\n\
[-reference-low-threshold | -ref-lt %d]    # values <= low threshold are not\n\
  considered\n\
[-reference-high-threshold | -ref-ht %d]   # values >= high threshold are not\n\
  considered\n\
[-reference-removed-fraction | -ref-rf %f] # maximal fraction of removed points\n\
  because of the threshold. If too many points are removed, the block is\n\
  discarded\n\
[-floating-selection-fraction[-ll|-hl] | -flo-frac[-ll|-hl] %lf] # fraction of\n\
  blocks from the floating image kept at a pyramid level, the blocks being\n\
  sorted w.r.t their variance (see note (1) for [-ll|-hl])\n\
### pairing ###\n\
[-search-neighborhood-half-size | -se-hsize %d %d %d] # half size of the search\n\
  neighborhood in the reference when looking fro similar blocks\n\
[-search-neighborhood-step | -se-step %d %d %d] # step between blocks to be\n\
  tested in the search neighborhood\n\
[-similarity-measure | -similarity | -si [cc|ecc|ssd|sad]]  # similarity measure\n\
  cc: correlation coefficient\n\
  ecc: extended correlation coefficient\n\
[-similarity-measure-threshold | -si-th %lf]    # threshold on the similarity\n\
  measure: pairings below that threshold are discarded\n			\
### transformation type ###\n\
[-transformation-type|-transformation|-trsf-type %s] # transformation type\n\
  translation2D, translation3D, translation-scaling2D, translation-scaling3D,\n\
  rigid2D, rigid3D, rigid, similitude2D, similitude3D, similitude,\n\
  affine2D, affine3D, affine, vectorfield2D, vectorfield3D, vectorfield, vector\n\
### transformation regularization ###\n\
[-elastic-regularization-sigma[-ll|-hl] | -elastic-sigma[-ll|-hl]  %lf %lf %lf]\n\
  # sigma for elastic regularization (only for vector field) (see note (1) for\n\
  [-ll|-hl])\n\
### transformation estimation ###\n\
[-estimator-type|-estimator|-es-type %s] # transformation estimator\n\
  wlts: weighted least trimmed squares\n\
  lts: least trimmed squares\n\
  wls: weighted least squares\n\
  ls: least squares\n\
[-lts-cut|-lts-fraction %lf] # for trimmed estimations, fraction of pairs that are kept\n\
[-lts-deviation %lf] # for trimmed estimations, defines the threshold to discard\n\
  pairings, ie 'average + this_value * standard_deviation'\n\
[-lts-iterations %d] # for trimmed estimations, the maximal number of iterations\n\
[-fluid-sigma|-lts-sigma[-ll|-hl] %lf %lf %lf] # sigma for fluid regularization,\n\
  ie field interpolation and regularization for pairings (only for vector field)\n\
  (see note (1) for [-ll|-hl])\n\
### end conditions for matching loop ###\n\
[-max-iteration[-ll|-hl]|-max-iter[-ll|-hl] %d]   # maximal number of iteration\n\
  (see note (1) for [-ll|-hl])\n\
[-corner-ending-condition|-rms] # evolution of image corners\n\
### filter type ###\n\
[-gaussian-filter-type|-filter-type deriche|fidrich|young-1995|young-2002|...\n\
  ...|gabor-young-2002|convolution] # type of filter for image/vector field\n\
  smoothing\n\
### parallelism ###\n\
[-parallel|-no-parallel] # use parallelism (or not)\n\
[-max-chunks %d] # maximal number of chunks for open mp\n\
[-parallel-scheduling|-ps default|static|dynamic-one|dynamic|guided] # type\n\
  of scheduling for open mp\n\
### misc writing stuff ###\n\
[-print-parameters] # print all parameters\n\
[-default-filenames|-df]     # use default filename names\n\
[-no-default-filenames|-ndf] # do not use default filename names\n\
[-command-line %s]           # write the command line\n\
[-logfile %s]                # write some output in this logfile\n\
[-vischeck]  # write an image with 'active' blocks\n\
[-write_def] # id. \n\
[-verbose] [-no-verbose]\n\
[-h|-help|--h|--help]\n\
\n\
Notes\n\
(1) If -ll or -hl are respectively added to the option, this specifies only the\n\
  value for respectively the lowest or the highest level of the pyramid (recall\n\
  that the most lowest level, ie #0, refers to the original image). For\n\
  intermediary levels, values are linearly interpolated.\n\
";





typedef struct local_parameter {
	
  /* file names
     - images
     - transformations
  */
  char *floating_image;
  char *reference_image;
  char *result_image;
  
  char *initial_real_transformation;
  char *initial_voxel_transformation;

  char *initial_result_real_transformation;
  char *initial_result_voxel_transformation;

  char *result_real_transformation;
  char *result_voxel_transformation;
 
  /* pre-processing before matching
   */
  int normalisation;


  /* parameters for hierarchical block matching
   */
  bal_blockmatching_pyramidal_param param;


  /* writing stuff
   */
  int use_default_filename;
  char *command_line_file;
  char *log_file;
  
  int print_time;

  /* misc
   */
  int print_parameters;


} local_parameter;



#define FILENAMELENGTH 128



static void _ErrorParse( char *str, int flag );
static void _Parse( int argc, char *argv[], local_parameter *p );
static void _InitParam( local_parameter *par );
static double _GetTime();
static double _GetClock();
static char *_BaseName( char *p );
#ifdef _UNUSED_
static void _WriteCommandLine( local_parameter *p, int argc, char *argv[] );
#endif






int main(int argc, char *argv[])
{
  int i;
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

  
  /* displaying help is required
   */
  i = 1;
  while ( i < argc ) {
    if ( ( strcmp ( argv[i], "-help") == 0 && argv[i][5] == '\0' ) 
	 || ( strcmp ( argv[i], "--help") == 0 && argv[i][6] == '\0' ) ) {
      _ErrorParse( NULL, 2 );
    }
    else if ( ( strcmp ( argv[i], "-h") == 0 && argv[i][2] == '\0' ) 
	      || ( strcmp ( argv[i], "--h") == 0 && argv[i][3] == '\0' ) ) {
      _ErrorParse( NULL, 1 );
    }
    i++;
  }

  
  /* parsing parameters 
   */
  _InitParam( &p );
  _Parse( argc, argv, &p );


   if( blockmatching(
		    p.floating_image,
		    p.reference_image,
		    p.result_image,
		    p.initial_real_transformation,
		    p.initial_voxel_transformation,
		    p.initial_result_real_transformation,
		    p.initial_result_voxel_transformation,
		    p.result_real_transformation,
		    p.result_voxel_transformation,
		    p.normalisation,
		    p.param,
		    p.use_default_filename,
		    p.command_line_file,
		    p.log_file,
		    p.print_time,
		    p.print_parameters,
		    _debug_
		     ) )
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

  return( 0 );
}







/************************************************************
 *
 *
 *
 ************************************************************/


static void _ErrorParse( char *str, int flag )
{
  FILE *output = stderr;

  if ( flag == 2 ) output = stdout;
  (void)fprintf( output, "Usage: %s %s", _BaseName(program), usage );
  if ( flag == 2 || flag == 1 )
    (void)fprintf( output, "%s\n", options );
  if ( flag == 2 )
    (void)fprintf( output, "%s\n", detail );
  if ( str != NULL ) (void)fprintf( output, "Error: %s\n", str );
  if ( flag != 2 && flag != 1 ) {
    (void)fprintf( output, "\n" );
    (void)fprintf( output, "# use '-h' or '--h' for a full list of options\n" );
    (void)fprintf( output, "# use '-help' or '--help' for more details\n" );
    (void)fprintf( output, "\n" );
  }
  exit( -1 );
}



/************************************************************
 *
 * reading parameters is done in two steps
 * 1. one looks for the transformation type
 *    -> this implies dedicated initialization
 * 2. the other parameters are read
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
    /* transformation type
     */
    if ( strcmp ( argv[i], "-transformation-type" ) == 0 
	 || strcmp ( argv[i], "-transformation" ) == 0
	 || strcmp ( argv[i], "-trsf-type" ) == 0 ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "-transformation-type", 0 );
      if ( strcmp ( argv[i], "translation2D" ) == 0 ) {
	p->param.transformation_type = TRANSLATION_2D;
      }
      else if ( strcmp ( argv[i], "translation3D" ) == 0 ) {
	p->param.transformation_type = TRANSLATION_3D;
      }
      else if ( strcmp ( argv[i], "translation" ) == 0 && argv[i][11] == '\0') {
	p->param.transformation_type = TRANSLATION_3D;
      }
      else if ( strcmp ( argv[i], "translation-scaling2D" ) == 0 ) {
	p->param.transformation_type = TRANSLATION_SCALING_2D;
      }
      else if ( strcmp ( argv[i], "translation-scaling3D" ) == 0 ) {
	p->param.transformation_type = TRANSLATION_SCALING_3D;
      }
      else if ( strcmp ( argv[i], "rigid2D" ) == 0 ) {
	p->param.transformation_type = RIGID_2D;
      }
      else if ( strcmp ( argv[i], "rigid3D" ) == 0 ) {
	p->param.transformation_type = RIGID_3D;
      }
      else if ( (strcmp ( argv[i], "rigid" ) == 0 && argv[i][5] == '\0') ) {
	p->param.transformation_type = RIGID_3D;
      }
      else if ( strcmp ( argv[i], "similitude2D" ) == 0 ) {
	p->param.transformation_type = SIMILITUDE_2D;
      }
      else if ( strcmp ( argv[i], "similitude3D" ) == 0 ) {
	p->param.transformation_type = SIMILITUDE_3D;
      }
      else if ( strcmp ( argv[i], "similitude" ) == 0 ) {
	p->param.transformation_type = SIMILITUDE_3D;
      }
      else if ( strcmp ( argv[i], "affine2D" ) == 0 ) {
	p->param.transformation_type = AFFINE_2D;
      }
      else if ( strcmp ( argv[i], "affine3D" ) == 0 ) {
	p->param.transformation_type = AFFINE_3D;
      }
      else if ( strcmp ( argv[i], "affine" ) == 0 ) {
	p->param.transformation_type = AFFINE_3D;
      }
      /*
	else if ( strcmp ( argv[i], "spline" ) == 0 ) {
	p->param.transformation_type = SPLINE;
	}
      */
      else if ( strcmp ( argv[i], "vectorfield" ) == 0 
		|| strcmp ( argv[i], "vector" ) == 0 ) {
	p->param.transformation_type = VECTORFIELD_3D;
      }
      else if ( strcmp ( argv[i], "vectorfield3D" ) == 0 
		|| strcmp ( argv[i], "vector3D" ) == 0 ) {
	p->param.transformation_type = VECTORFIELD_3D;
      }
      else if ( strcmp ( argv[i], "vectorfield2D" ) == 0 
		|| strcmp ( argv[i], "vector2D" ) == 0 ) {
	p->param.transformation_type = VECTORFIELD_2D;
      }
      else {
	fprintf( stderr, "unknown transformation type: '%s'\n", argv[i] );
	_ErrorParse( "-transformation-type", 0 );
      }
    }
  }



  /* initialisation of parameters dedicated to either
     linear transformation or vector field
  */

  switch ( p->param.transformation_type ) {
  default :
    _ErrorParse( "unknown transformation type", 0 );
    break;
  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :
    BAL_InitBlockMatchingPyramidalParametersForLinearTransformation( &(p->param) );
    break;
  case VECTORFIELD_3D :
  case VECTORFIELD_2D :
    BAL_InitBlockMatchingPyramidalParametersForVectorfieldTransformation( &(p->param) );
    break;
  }
  


  /* reading the other parameters
   */

  for ( i=1; i<argc; i++ ) {
		
    /* transformation type
       already read
    */
    if ( strcmp ( argv[i], "-transformation-type" ) == 0 
	 || strcmp ( argv[i], "-transformation" ) == 0
	 || strcmp ( argv[i], "-trsf-type" ) == 0 ) {
      i ++;
    }
    
    /* image file names 
     */
    else if ( strcmp ( argv[i], "-reference") == 0
	      || (strcmp ( argv[i], "-ref") == 0 && argv[i][4] == '\0') ) {
      i++;
      if ( i >= argc) _ErrorParse( "parsing -reference", 0 );
      p->reference_image = argv[i];
    }
    else if ( strcmp ( argv[i], "-floating") == 0
	      || (strcmp ( argv[i], "-flo") == 0 && argv[i][4] == '\0') ) {
      i++;
      if ( i >= argc) _ErrorParse( "parsing -floating", 0 );
      p->floating_image = argv[i];
    }
    else if ( strcmp ( argv[i], "-result") == 0
	      || (strcmp ( argv[i], "-res") == 0 && argv[i][4] == '\0') ) {
      i++;
      if ( i >= argc) _ErrorParse( "parsing -result", 0 );
      p->result_image = argv[i];
    }

    /* transformation file names 
     */
    else if ( strcmp ( argv[i], "-initial-transformation" ) == 0
	      || strcmp ( argv[i], "-init-trsf" ) == 0 ) {
      i++;
      if ( i >= argc) _ErrorParse( "parsing -initial-transformation", 0 );
      p->initial_real_transformation = argv[i];
    }
    else if ( strcmp ( argv[i], "-initial-voxel-transformation" ) == 0
	      || strcmp ( argv[i], "-init-voxel-trsf" ) == 0 ) {
      i++;
      if ( i >= argc) _ErrorParse( "parsing -initial-voxel-transformation", 0 );
      p->initial_voxel_transformation = argv[i];
    }

    else if ( strcmp ( argv[i], "-initial-result-transformation" ) == 0
	      || strcmp ( argv[i], "-init-res-trsf" ) == 0 ) {
      i++;
      if ( i >= argc) _ErrorParse( "parsing -initial-result-transformation", 0 );
      p->initial_result_real_transformation = argv[i];
    }
    else if ( strcmp ( argv[i], "-initial-result-voxel-transformation" ) == 0
	      || strcmp ( argv[i], "-init-res-voxel-trsf" ) == 0 ) {
      i++;
      if ( i >= argc) _ErrorParse( "parsing -initial-result-voxel-transformation", 0 );
      p->initial_result_voxel_transformation = argv[i];
    }

    else if ( strcmp ( argv[i], "-result-transformation" ) == 0
	      || strcmp ( argv[i], "-res-trsf" ) == 0 ) {
      i++;
      if ( i >= argc) _ErrorParse( "parsing -result-transformation", 0 );
      p->result_real_transformation = argv[i];
    }
    else if ( strcmp ( argv[i], "-result-voxel-transformation" ) == 0
	      || strcmp ( argv[i], "-res-voxel-trsf" ) == 0 ) {
      i++;
      if ( i >= argc) _ErrorParse( "parsing -result-voxel-transformation", 0 );
      p->result_voxel_transformation = argv[i];
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
      status = sscanf( argv[i], "%d", &(p->param.pyramid_lowest_level) );
      if ( status <= 0 ) _ErrorParse( "parsing -pyramid-lowest-level", 0 );
    }
    else if ( strcmp ( argv[i], "-pyramid-highest-level" ) == 0
	      || (strcmp( argv[i], "-py-hl") == 0 && argv[i][6] == '\0') ) { 
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -pyramid-highest-level", 0 );
      status = sscanf( argv[i], "%d", &(p->param.pyramid_highest_level) );
      if ( status <= 0 ) _ErrorParse( "parsing -pyramid-highest-level", 0 );
    }
    else if ( strcmp ( argv[i], "-pyramid-gaussian-filtering" ) == 0
	      || (strcmp( argv[i], "-py-gf") == 0 && argv[i][6] == '\0') ) { 
      p->param.pyramid_gaussian_filtering = 1;
    }
    

    /* block geometry
     */
    else if ( strcmp (argv[i], "-block-size" ) == 0 
	      || strcmp (argv[i], "-bl-size" ) == 0 ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -block-size %d", 0 );
      status = sscanf( argv[i], "%d", &(p->param.block_dim.x) );
      if ( status <= 0 ) _ErrorParse( "parsing -block-size %d", 0 );
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -block-size %d %d", 0 );
      status = sscanf( argv[i], "%d", &(p->param.block_dim.y) );
      if ( status <= 0 ) _ErrorParse( "parsing -block-sizes %d %d", 0 );
      i ++;
      if ( i >= argc) p->param.block_dim.z = 1;
      else {
	status = sscanf( argv[i], "%d", &(p->param.block_dim.z) );
	if ( status <= 0 ) {
	  i--;
	  p->param.block_dim.z = 1;
	}
      }
    }
    else if ( strcmp (argv[i], "-block-spacing" ) == 0 
	      || strcmp (argv[i], "-bl-space") == 0 ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -block-spacing %d", 0 );
      status = sscanf( argv[i], "%d", &(p->param.block_spacing.x) );
      if ( status <= 0 ) _ErrorParse( "parsing -block-spacing %d", 0 );
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -block-spacing %d %d", 0 );
      status = sscanf( argv[i], "%d", &(p->param.block_spacing.y) );
      if ( status <= 0 ) _ErrorParse( "parsing -block-spacing %d %d", 0 );
      i ++;
      if ( i >= argc) p->param.block_spacing.z = 0;
      else {
	status = sscanf( argv[i], "%d", &(p->param.block_spacing.z) );
	if ( status <= 0 ) {
	  i--;
	  p->param.block_spacing.z = 0;
	}
      }
    }
    else if ( strcmp (argv[i], "-block-border" ) == 0 
	      || strcmp (argv[i], "-bl-border") == 0  ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -block-border %d", 0 );
      status = sscanf( argv[i], "%d", &(p->param.block_border.x) );
      if ( status <= 0 ) _ErrorParse( "parsing -block-borders %d", 0 );
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -block-borders %d %d", 0 );
      status = sscanf( argv[i], "%d", &(p->param.block_border.y) );
      if ( status <= 0 ) _ErrorParse( "parsing -block-borders %d %d", 0 );
      i ++;
      if ( i >= argc) p->param.block_border.z = 0;
      else {
	status = sscanf( argv[i], "%d", &(p->param.block_border.z) );
	if ( status <= 0 ) {
	  i--;
	  p->param.block_border.z = 0;
	}
      }
    }
    

    /* block selection
     */
    else if ( strcmp ( argv[i], "-floating-low-threshold") == 0 
	      || strcmp ( argv[i], "-flo-lt") == 0 ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -floating-low-threshold", 0 );
      status = sscanf( argv[i], "%d", &(p->param.floating_selection.low_threshold) );
      if ( status <= 0 ) _ErrorParse( "parsing -floating-low-threshold", 0 );
    }
    else if ( strcmp ( argv[i], "-floating-high-threshold") == 0 
	      || strcmp ( argv[i], "-flo-ht") == 0 ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -floating-high-threshold", 0 );
      status = sscanf( argv[i], "%d", &(p->param.floating_selection.high_threshold) );
      if ( status <= 0 ) _ErrorParse( "parsing -floating-high-threshold", 0 );
    }
    else if ( strcmp ( argv[i], "-floating-removed-fraction") == 0 
	      || strcmp ( argv[i], "-flo-rf") == 0 ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -floating-removed-fraction", 0 );
      status = sscanf( argv[i], "%lf", &(p->param.floating_selection.max_removed_fraction) );
      if ( status <= 0 ) _ErrorParse( "parsing -floating-removed-fraction", 0 );
    }

    else if ( strcmp ( argv[i], "-reference-low-threshold") == 0 
	      || strcmp ( argv[i], "-ref-lt") == 0 ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -reference-low-threshold", 0 );
      status = sscanf( argv[i], "%d", &(p->param.reference_selection.low_threshold) );
      if ( status <= 0 ) _ErrorParse( "parsing -reference-low-threshold", 0 );
    }
    else if ( strcmp ( argv[i], "-reference-high-threshold") == 0 
	      || strcmp ( argv[i], "-ref-ht") == 0 ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -reference-high-threshold", 0 );
      status = sscanf( argv[i], "%d", &(p->param.reference_selection.high_threshold) );
      if ( status <= 0 ) _ErrorParse( "parsing -reference-high-threshold", 0 );
    }
    else if ( strcmp ( argv[i], "-reference-removed-fraction") == 0 
	      || strcmp ( argv[i], "-ref-rf") == 0 ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -reference-removed-fraction", 0 );
      status = sscanf( argv[i], "%lf", &(p->param.reference_selection.max_removed_fraction) );
      if ( status <= 0 ) _ErrorParse( "parsing -reference-removed-fraction", 0 );
    }

    else if ( ( strcmp ( argv[i], "-floating-selection-fraction") == 0 && argv[i][28] == '\0' )
	      || ( strcmp ( argv[i], "-flo-frac") == 0 && argv[i][9] == '\0' ) ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "-floating-selection-fraction", 0 );
      status = sscanf( argv[i], "%lf", &(p->param.blocks_fraction.lowest) );
      if ( status <= 0 ) _ErrorParse( "-floating-selection-fraction", 0 );
      p->param.blocks_fraction.highest = p->param.blocks_fraction.lowest;
    }
    else if ( strcmp ( argv[i], "-floating-selection-fraction-lowest-level") == 0 
	      || strcmp ( argv[i], "-floating-selection-fraction-ll") == 0 
	      || strcmp ( argv[i], "-flo-frac-ll") == 0 ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "-floating-selection-fraction-lowest-level", 0 );
      status = sscanf( argv[i], "%lf", &(p->param.blocks_fraction.lowest) );
      if ( status <= 0 ) _ErrorParse( "-floating-selection-fraction-lowest-level", 0 );
    }
    else if ( strcmp ( argv[i], "-floating-selection-fraction-highest-level") == 0 
	      || strcmp ( argv[i], "-floating-selection-fraction-hl") == 0 
	      || strcmp ( argv[i], "-flo-frac-hl") == 0 ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "-floating-selection-fraction-highest-level", 0 );
      status = sscanf( argv[i], "%lf", &(p->param.blocks_fraction.highest) );
      if ( status <= 0 ) _ErrorParse( "-floating-selection-fraction-highest-level", 0 );
    }


    /* pairing
     */
    else if ( strcmp (argv[i], "-search-neighborhood-half-size" ) == 0 
	      || strcmp (argv[i], "-se-hsize" ) == 0 ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -search-neighborhood-half-size %d", 0 );
      status = sscanf( argv[i], "%d", &(p->param.half_neighborhood_size.x) );
      if ( status <= 0 ) _ErrorParse( "parsing -search-neighborhood-half-size %d", 0 );
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -search-neighborhood-half-size %d %d", 0 );
      status = sscanf( argv[i], "%d", &(p->param.half_neighborhood_size.y) );
      if ( status <= 0 ) _ErrorParse( "parsing -search-neighborhood-half-size %d %d", 0 );
      i ++;
      if ( i >= argc) p->param.half_neighborhood_size.z = 0;
      else {
	status = sscanf( argv[i], "%d", &(p->param.half_neighborhood_size.z) );
	if ( status <= 0 ) {
	  i--;
	  p->param.half_neighborhood_size.z = 0;
	}
      }
    }
    else if ( strcmp (argv[i], "-search-neighborhood-step" ) == 0 
	      || strcmp (argv[i], "-se-step") == 0 ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -search-neighborhood-step %d", 0 );
      status = sscanf( argv[i], "%d", &(p->param.step_neighborhood_search.x) );
      if ( status <= 0 ) _ErrorParse( "parsing -search-neighborhood-step %d", 0 );
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -search-neighborhood-step %d %d", 0 );
      status = sscanf( argv[i], "%d", &(p->param.step_neighborhood_search.y) );
      if ( status <= 0 ) _ErrorParse( "parsing -search-neighborhood-step %d %d", 0 );
      i ++;
      if ( i >= argc) p->param.step_neighborhood_search.z = 0;
      else {
	status = sscanf( argv[i], "%d", &(p->param.step_neighborhood_search.z) );
	if ( status <= 0 ) {
	  i--;
	  p->param.step_neighborhood_search.z = 0;
	}
      }
    }

    else if ( strcmp ( argv[i], "-similarity-measure" ) == 0
	      || strcmp ( argv[i], "-similarity" ) == 0 
	      || (strcmp ( argv[i], "-si" ) == 0 && argv[i][3] == '\0') ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "-similarity-measure", 0 );
      if ( strcmp ( argv[i], "ssd" ) == 0 && argv[i][3] == '\0' ) {
	p->param.similarity_measure = _SSD_;
      }
      else if ( strcmp ( argv[i], "sad" ) == 0 && argv[i][3] == '\0' ) {
	p->param.similarity_measure = _SAD_;
      }
      else if ( strcmp ( argv[i], "cc" ) == 0 && argv[i][2] == '\0' ) {
	p->param.similarity_measure = _SQUARED_CC_;
      }
      else if ( strcmp ( argv[i], "ecc" ) == 0 && argv[i][3] == '\0' ) {
	p->param.similarity_measure = _SQUARED_EXTCC_;
      }
      else {
	fprintf( stderr, "unknown similarity measure: '%s'\n", argv[i] );
	_ErrorParse( "-similarity-measure", 0 );
      }
    }
    
    else if ( strcmp ( argv[i], "-similarity-measure-threshold" ) == 0
	      || strcmp ( argv[i], "-si-th" ) == 0 ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "-similarity-measure-threshold", 0 );
      status = sscanf( argv[i], "%lf", &(p->param.similarity_measure_threshold) );
      if ( status <= 0 ) _ErrorParse( "-similarity-measure-threshold", 0 );
    }
    


    /* transformation definition and computation 
     */
    else if ( ( strcmp ( argv[i], "-elastic-regularization-sigma" ) == 0  && argv[i][29] == '\0' )
	      || ( strcmp (argv[i], "-elastic-sigma" ) == 0 && argv[i][14] == '\0' ) ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -elastic-regularization-sigma %lf", 0 );
      status = sscanf( argv[i], "%lf", &(p->param.elastic_regularization_sigma.lowest.x) );
      if ( status <= 0 ) _ErrorParse( "parsing -elastic-regularization-sigma %lf", 0 );
      i ++;
      if ( i >= argc) {
	p->param.elastic_regularization_sigma.lowest.y = p->param.elastic_regularization_sigma.lowest.x;
	p->param.elastic_regularization_sigma.lowest.z = p->param.elastic_regularization_sigma.lowest.x;
      }
      else {
	status = sscanf( argv[i], "%lf", &(p->param.elastic_regularization_sigma.lowest.y) );
	if ( status <= 0 ) {
	  i--;
	  p->param.elastic_regularization_sigma.lowest.y = p->param.elastic_regularization_sigma.lowest.x;
	  p->param.elastic_regularization_sigma.lowest.z = p->param.elastic_regularization_sigma.lowest.x;
	}
	else {
	  i ++;
	  if ( i >= argc) p->param.elastic_regularization_sigma.lowest.z = 0;
	  else {
	    status = sscanf( argv[i], "%lf", &(p->param.elastic_regularization_sigma.lowest.z) );
	    if ( status <= 0 ) {
	      i--;
	      p->param.elastic_regularization_sigma.lowest.z = 0;
	    }
	  }
	}
      }
      p->param.elastic_regularization_sigma.highest.x = p->param.elastic_regularization_sigma.lowest.x;
      p->param.elastic_regularization_sigma.highest.y = p->param.elastic_regularization_sigma.lowest.y;
      p->param.elastic_regularization_sigma.highest.z = p->param.elastic_regularization_sigma.lowest.z;
    }
    else if ( strcmp ( argv[i], "-elastic-regularization-sigma-lowest-level" ) == 0  
	      || strcmp ( argv[i], "-elastic-regularization-sigma-ll" ) == 0  
	      || strcmp (argv[i], "-elastic-sigma-ll" ) == 0 ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -elastic-regularization-sigma-lowest-level %lf", 0 );
      status = sscanf( argv[i], "%lf", &(p->param.elastic_regularization_sigma.lowest.x) );
      if ( status <= 0 ) _ErrorParse( "parsing -elastic-regularization-sigma-lowest-level %lf", 0 );
      i ++;
      if ( i >= argc) {
	p->param.elastic_regularization_sigma.lowest.y = p->param.elastic_regularization_sigma.lowest.x;
	p->param.elastic_regularization_sigma.lowest.z = p->param.elastic_regularization_sigma.lowest.x;
      }
      else {
	status = sscanf( argv[i], "%lf", &(p->param.elastic_regularization_sigma.lowest.y) );
	if ( status <= 0 ) {
	  i--;
	  p->param.elastic_regularization_sigma.lowest.y = p->param.elastic_regularization_sigma.lowest.x;
	  p->param.elastic_regularization_sigma.lowest.z = p->param.elastic_regularization_sigma.lowest.x;
	}
	else {
	  i ++;
	  if ( i >= argc) p->param.elastic_regularization_sigma.lowest.z = 0;
	  else {
	    status = sscanf( argv[i], "%lf", &(p->param.elastic_regularization_sigma.lowest.z) );
	    if ( status <= 0 ) {
	      i--;
	      p->param.elastic_regularization_sigma.lowest.z = 0;
	    }
	  }
	}
      }
    }
    else if ( strcmp ( argv[i], "-elastic-regularization-sigma-highest-level" ) == 0  
	      || strcmp ( argv[i], "-elastic-regularization-sigma-hl" ) == 0  
	      || strcmp (argv[i], "-elastic-sigma-hl" ) == 0 ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -elastic-regularization-sigma-highest-level %lf", 0 );
      status = sscanf( argv[i], "%lf", &(p->param.elastic_regularization_sigma.highest.x) );
      if ( status <= 0 ) _ErrorParse( "parsing -elastic-regularization-sigma-highest-level %lf", 0 );
      i ++;
      if ( i >= argc) {
	p->param.elastic_regularization_sigma.highest.y = p->param.elastic_regularization_sigma.highest.x;
	p->param.elastic_regularization_sigma.highest.z = p->param.elastic_regularization_sigma.highest.x;
      }
      else {
	status = sscanf( argv[i], "%lf", &(p->param.elastic_regularization_sigma.highest.y) );
	if ( status <= 0 ) {
	  i--;
	  p->param.elastic_regularization_sigma.highest.y = p->param.elastic_regularization_sigma.highest.x;
	  p->param.elastic_regularization_sigma.highest.z = p->param.elastic_regularization_sigma.highest.x;
	}
	else {
	  i ++;
	  if ( i >= argc) p->param.elastic_regularization_sigma.highest.z = 0;
	  else {
	    status = sscanf( argv[i], "%lf", &(p->param.elastic_regularization_sigma.highest.z) );
	    if ( status <= 0 ) {
	      i--;
	      p->param.elastic_regularization_sigma.highest.z = 0;
	    }
	  }
	}
      }
    }



    /* estimator definition and computation 
     */    		
    else if ( strcmp ( argv[i], "-estimator-type") == 0 
	      || strcmp ( argv[i], "-estimator") == 0
	      || strcmp ( argv[i], "-es-type") == 0 ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "-estimator-type", 0 );
      if ( (strcmp ( argv[i], "ltsw" ) == 0 && argv[i][4] == '\0')
	   || (strcmp ( argv[i], "wlts" ) == 0 && argv[i][4] == '\0') ) {
	p->param.estimator.lowest.type = TYPE_WLTS;
      }
      else if ( strcmp ( argv[i], "lts" ) == 0 && argv[i][3] == '\0' ) {
	p->param.estimator.lowest.type = TYPE_LTS;
      }
      else if ( (strcmp ( argv[i], "lsw" ) == 0 && argv[i][3] == '\0')
		|| (strcmp ( argv[i], "wls" ) == 0 && argv[i][3] == '\0') ) {
	p->param.estimator.lowest.type = TYPE_WLS;
      }
      else if ( strcmp ( argv[i], "ls" ) == 0 && argv[i][2] == '\0' ) {
	p->param.estimator.lowest.type = TYPE_LS;
      }
      else {
	fprintf( stderr, "unknown estimator type: '%s'\n", argv[i] );
	_ErrorParse( "-estimator-type", 0 );
      }
      p->param.estimator.highest.type = p->param.estimator.lowest.type;
    }
		
    else if ( strcmp ( argv[i], "-lts-fraction" ) == 0 
	      || strcmp ( argv[i], "-lts-cut" ) == 0) {
      i ++;
      if ( i >= argc)    _ErrorParse( "-lts-fraction", 0 );
      status = sscanf( argv[i], "%lf", &(p->param.estimator.lowest.retained_fraction) );
      if ( status <= 0 ) _ErrorParse( "-lts-fraction", 0 );
      p->param.estimator.highest.retained_fraction = p->param.estimator.lowest.retained_fraction;
    }

    else if ( strcmp ( argv[i], "-lts-deviation" ) == 0 ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "-lts-deviation", 0 );
      status = sscanf( argv[i], "%lf", &(p->param.estimator.lowest.standard_deviation_threshold) );
      if ( status <= 0 ) _ErrorParse( "-lts-deviation", 0 );
      p->param.estimator.highest.standard_deviation_threshold = p->param.estimator.lowest.standard_deviation_threshold;
    }

    else if ( strcmp ( argv[i], "-lts-iterations" ) == 0 ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "-lts-iterations", 0 );
      status = sscanf( argv[i], "%d", &(p->param.estimator.lowest.max_iterations) );
      if ( status <= 0 ) _ErrorParse( "-lts-iterations", 0 );
      p->param.estimator.highest.max_iterations = p->param.estimator.lowest.max_iterations;
    }
    
    else if ( strcmp (argv[i], "-fluid-sigma-lowest-level" ) == 0
	      || strcmp (argv[i], "-fluid-sigma-ll" ) == 0
	      || strcmp (argv[i], "-lts-sigma-lowest-level" ) == 0
	      || strcmp (argv[i], "-lts-sigma-ll" ) == 0  ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -lts-sigma-lowest-level %lf", 0 );
      status = sscanf( argv[i], "%lf", &(p->param.estimator.lowest.sigma.x) );
      if ( status <= 0 ) _ErrorParse( "parsing -lts-sigma-lowest-level %lf", 0 );
      i ++;
      if ( i >= argc) {
	p->param.estimator.lowest.sigma.y = p->param.estimator.lowest.sigma.x;
	p->param.estimator.lowest.sigma.z = p->param.estimator.lowest.sigma.x;
      }
      else {
	status = sscanf( argv[i], "%lf", &(p->param.estimator.lowest.sigma.y) );
	if ( status <= 0 ) {
	  i--;
	  p->param.estimator.lowest.sigma.y = p->param.estimator.lowest.sigma.x;
	  p->param.estimator.lowest.sigma.z = p->param.estimator.lowest.sigma.x;
	}
	else {
	  i ++;
	  if ( i >= argc) p->param.estimator.lowest.sigma.z = 0;
	  else {
	    status = sscanf( argv[i], "%lf", &(p->param.estimator.lowest.sigma.z) );
	    if ( status <= 0 ) {
	      i--;
	      p->param.estimator.lowest.sigma.z = 0;
	    }
	  }
	}
      }
    }
    else if ( strcmp (argv[i], "-fluid-sigma-highest-level" ) == 0
	      || strcmp (argv[i], "-fluid-sigma-hl" ) == 0  
	      || strcmp (argv[i], "-lts-sigma-highest-level" ) == 0
	      || strcmp (argv[i], "-lts-sigma-hl" ) == 0) {
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -lts-sigma-highest-level %lf", 0 );
      status = sscanf( argv[i], "%lf", &(p->param.estimator.highest.sigma.x) );
      if ( status <= 0 ) _ErrorParse( "parsing -lts-sigma-highest-level %lf", 0 );
      i ++;
      if ( i >= argc) {
	p->param.estimator.highest.sigma.y = p->param.estimator.highest.sigma.x;
	p->param.estimator.highest.sigma.z = p->param.estimator.highest.sigma.x;
      }
      else {
	status = sscanf( argv[i], "%lf", &(p->param.estimator.highest.sigma.y) );
	if ( status <= 0 ) {
	  i--;
	  p->param.estimator.highest.sigma.y = p->param.estimator.highest.sigma.x;
	  p->param.estimator.highest.sigma.z = p->param.estimator.highest.sigma.x;
	}
	else {
	  i ++;
	  if ( i >= argc) p->param.estimator.highest.sigma.z = 0;
	  else {
	    status = sscanf( argv[i], "%lf", &(p->param.estimator.highest.sigma.z) );
	    if ( status <= 0 ) {
	      i--;
	      p->param.estimator.highest.sigma.z = 0;
	    }
	  }
	}
      }
    }
    else if ( (strcmp (argv[i], "-fluid-sigma" ) == 0 && argv[i][12] == '\0')
	      || (strcmp (argv[i], "-lts-sigma" ) == 0 && argv[i][10] == '\0') ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "parsing -lts-sigma %lf", 0 );
      status = sscanf( argv[i], "%lf", &(p->param.estimator.lowest.sigma.x) );
      if ( status <= 0 ) _ErrorParse( "parsing -lts-sigma %lf", 0 );
      i ++;
      if ( i >= argc) {
	p->param.estimator.lowest.sigma.y = p->param.estimator.lowest.sigma.x;
	p->param.estimator.lowest.sigma.z = p->param.estimator.lowest.sigma.x;
      }
      else {
	status = sscanf( argv[i], "%lf", &(p->param.estimator.lowest.sigma.y) );
	if ( status <= 0 ) {
	  i--;
	  p->param.estimator.lowest.sigma.y = p->param.estimator.lowest.sigma.x;
	  p->param.estimator.lowest.sigma.z = p->param.estimator.lowest.sigma.x;
	}
	else {
	  i ++;
	  if ( i >= argc) p->param.estimator.lowest.sigma.z = 0;
	  else {
	    status = sscanf( argv[i], "%lf", &(p->param.estimator.lowest.sigma.z) );
	    if ( status <= 0 ) {
	      i--;
	      p->param.estimator.lowest.sigma.z = 0;
	    }
	  }
	}
      }
      p->param.estimator.highest.sigma = p->param.estimator.lowest.sigma;
    }



    /* end conditions for matching loop
     */
    else if ( strcmp ( argv[i], "-max-iteration-lowest-level" ) == 0
	      || strcmp ( argv[i], "-max-iteration-ll" ) == 0
	      || strcmp ( argv[i], "-max-iterations-lowest-level" ) == 0
	      || strcmp ( argv[i], "-max-iterations-ll" ) == 0
	      || strcmp( argv[i], "-max-iter-lowest-level" ) == 0  
	      || strcmp( argv[i], "-max-iter-ll" ) == 0  ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "-max-iteration-lowest-level", 0 );
      status = sscanf( argv[i], "%d", &(p->param.max_iterations.lowest) );
      if ( status <= 0 ) _ErrorParse( "-max-iteration-lowest-level", 0 );
    }
    else if ( strcmp ( argv[i], "-max-iteration-highest-level" ) == 0
	      || strcmp ( argv[i], "-max-iteration-hl" ) == 0
	      || strcmp ( argv[i], "-max-iterations-highest-level" ) == 0
	      || strcmp ( argv[i], "-max-iterations-hl" ) == 0
	      || strcmp( argv[i], "-max-iter-highest-level" ) == 0  
	      || strcmp( argv[i], "-max-iter-hl" ) == 0  ) {
      i ++;
      if ( i >= argc)    _ErrorParse( "-max-iteration-highest-level", 0 );
      status = sscanf( argv[i], "%d", &(p->param.max_iterations.highest) );
      if ( status <= 0 ) _ErrorParse( "-max-iteration-highest-level", 0 );
    }
    else if ( (strcmp ( argv[i], "-max-iteration" ) == 0 && argv[i][14] == '\0')
	      || (strcmp ( argv[i], "-max-iterations" ) == 0 && argv[i][15] == '\0')
	      || (strcmp( argv[i], "-max-iter" ) == 0 && argv[i][9] == '\0') ) {
      i ++;                         
      if ( i >= argc)    _ErrorParse( "-max-iteration", 0 );
      status = sscanf( argv[i], "%d", &(p->param.max_iterations.lowest) );
      if ( status <= 0 ) _ErrorParse( "-max-iteration", 0 );
      p->param.max_iterations.highest = p->param.max_iterations.lowest;
    }
		
    else if ( strcmp ( argv[i], "-corner-ending-condition" ) == 0
	      || (strcmp ( argv[i], "-rms") == 0 && argv[i][4] == '\0') ) {
      p->param.rms_ending_condition = 1;
    }

    
    /* filter type for image smoothing
     */
    else if ( strcmp ( argv[i], "-gaussian-filter-type" ) == 0 
	      || strcmp ( argv[i], "-filter-type" ) == 0 ) {
      i++;
      if ( i >= argc)    _ErrorParse( "-gaussian-filter-type", 0 );
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
    else if ( strcmp ( argv[i], "-print-parameters" ) == 0 ) {
      p->print_parameters = 1;
    }				

    else if ( strcmp ( argv[i], "-default-filenames" ) == 0 
	      || (strcmp ( argv[i], "-df" ) == 0 && argv[i][3] == '\0') ) {
      p->use_default_filename = 1;
    }
    else if ( strcmp ( argv[i], "-no-default-filenames" ) == 0 
	      || (strcmp ( argv[i], "-ndf" ) == 0 && argv[i][4] == '\0') ) {
      p->use_default_filename = 0;
    }

    else if ( strcmp ( argv[i], "-command-line") == 0 ) {
      i++;
      if ( i >= argc) _ErrorParse( "parsing -command-line", 0 );
      p->command_line_file = argv[i];
    }
    else if ( strcmp ( argv[i], "-logfile" ) == 0  ) {
      i++;
      if ( i >= argc) _ErrorParse( "parsing -logfile", 0 );
      p->log_file = argv[i];
      if ( p->param.verbose <= 0 ) p->param.verbose = 1;
    }

    else if ( (strcmp ( argv[i], "-time" ) == 0 && argv[i][5] == '\0') ) {
      p->print_time = 1;
    }
    else if ( (strcmp ( argv[i], "-notime" ) == 0 && argv[i][7] == '\0')  
	      || (strcmp ( argv[i], "-no-time" ) == 0 && argv[i][8] == '\0') ) {
      p->print_time = 0;
    }


    else if (strcmp ( argv[i], "-vischeck") == 0){
      p->param.vischeck = 1;
    }
    else if (strcmp ( argv[i], "-write_def") == 0){
      p->param.write_def = 1;
    }
    
    else if ( (strcmp ( argv[i], "-verbose") == 0 && argv[i][8] == '\0')
	      || (strcmp ( argv[i], "-v") == 0 && argv[i][2] == '\0') ) {
      if ( p->param.verbose <= 0 ) p->param.verbose = 1;
      else p->param.verbose ++;

      BAL_IncrementVerboseInBalBlockTools(  );
      BAL_IncrementVerboseInBalBlock(  );
      BAL_IncrementVerboseInBalEstimator(  );
      BAL_IncrementVerboseInBalFieldTools(  );
      BAL_IncrementVerboseInBalImageTools(  );
      BAL_IncrementVerboseInBalImage(  );
      BAL_IncrementVerboseInBalLinearTrsf(  );
      BAL_IncrementVerboseInBalMatching(  );
      BAL_IncrementVerboseInBalMatrix(  );
      BAL_IncrementVerboseInBalPyramid(  );
      BAL_IncrementVerboseInBalTransformationTools(  );
      BAL_IncrementVerboseInBalTransformation(  );
      BAL_IncrementVerboseInBalVectorField(  );

      incrementVerboseInChunks(  );
      incrementVerboseInReech4x4();
      incrementVerboseInReechDef();
      Recline_verbose();

    }
    else if ( (strcmp ( argv[i], "-no-verbose") == 0) 
	      || (strcmp ( argv[i], "-nv") == 0 && argv[i][3] == '\0') ) {
      p->param.verbose = 0;
      
      BAL_DecrementVerboseInBalBlockTools(  );
      BAL_DecrementVerboseInBalBlock(  );
      BAL_DecrementVerboseInBalEstimator(  );
      BAL_DecrementVerboseInBalFieldTools(  );
      BAL_DecrementVerboseInBalImageTools(  );
      BAL_DecrementVerboseInBalImage(  );
      BAL_DecrementVerboseInBalLinearTrsf(  );
      BAL_DecrementVerboseInBalMatching(  );
      BAL_DecrementVerboseInBalMatrix(  );
      BAL_DecrementVerboseInBalPyramid(  );
      BAL_DecrementVerboseInBalTransformationTools(  );
      BAL_DecrementVerboseInBalTransformation(  );
      BAL_DecrementVerboseInBalVectorField(  );

      decrementVerboseInChunks(  );
      decrementVerboseInReech4x4();
      decrementVerboseInReechDef();
      Recline_noverbose();

    }

    else if ( (strcmp ( argv[i], "-debug") == 0 && argv[i][6] == '\0')
	      || (strcmp ( argv[i], "-D") == 0 && argv[i][2] == '\0') ) {

      BAL_IncrementDebugInBalFieldTools(  );
      BAL_IncrementDebugInBalMatching(  );
      incrementDebugInChunks(  );
    }

    else if ( strcmp ( argv[i], "-no-debug") == 0 && argv[i][9] == '\0' ) {

      BAL_DecrementDebugInBalFieldTools(  );
      BAL_DecrementDebugInBalMatching(  );
      decrementDebugInChunks(  );

    }

    else if ( strcmp ( argv[i], "-trace") == 0 && argv[i][6] == '\0' ) {

      BAL_IncrementTraceInBalMatching(  );

    }

    else if ( strcmp ( argv[i], "-no-trace") == 0 && argv[i][9] == '\0' ) {

      BAL_DecrementTraceInBalMatching(  );

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
}





static void _InitParam( local_parameter *p ) 
{
  /* file names
     - images
     - transformations
  */
  p->floating_image = NULL;
  p->reference_image = NULL;
  p->result_image = NULL;
  
  p->initial_real_transformation = NULL;
  p->initial_voxel_transformation = NULL;

  p->initial_result_real_transformation = NULL;
  p->initial_result_voxel_transformation = NULL;

  p->result_real_transformation = NULL;
  p->result_voxel_transformation = NULL;
  

  /* pre-processing before matching
   */
#ifdef _ORIGINAL_BALADIN_UNSIGNED_CHAR_IMAGE_
  p->normalisation = 1;
#else
  p->normalisation = 0;
#endif


  /* parameters for hierarchical block matching
   */
  BAL_InitBlockMatchingPyramidalParameters( &(p->param) );


  /* writing stuff
   */
  p->use_default_filename = 1;
  p->command_line_file = NULL;
  p->log_file = NULL;
  p->print_time = 1;
  
  /* misc
   */
  p->print_parameters = 0;

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



#ifdef _UNUSED_
static void _WriteCommandLine( local_parameter *p, int argc, char *argv[] )
{
  char *proc="_WriteCommandLine";
  FILE *f = NULL;
  int i, j;
  char filename[FILENAMELENGTH];
  time_t t = time(NULL);

  if ( p->command_line_file != NULL ) { 
    f = fopen( p->command_line_file, "w" );
    if ( f == NULL ) 
      fprintf( stderr, "%s: unable to open '%s' for writing\n", proc, p->command_line_file );
  }
  else if ( p->use_default_filename == 1 ) {
    sprintf( filename, "%s-%d.sh", _BaseName( argv[0] ), getpid() );
    f = fopen( filename, "w" );
    if ( f == NULL ) 
      fprintf( stderr, " unable to open '%s' for writing\n", filename );
  }
  
  if ( f != NULL ) {
    fprintf( f, "%% %s\n", ctime( &t ) ); 
    for (i=0; i<argc; i++) {
      for ( j=0; j<strlen( argv[i] ); j++ )
	fprintf( f, "%c", argv[i][j] );
      fprintf( f, " " );
    }
    fprintf( f, "\n" );
    fclose( f );
  }

}
#endif
