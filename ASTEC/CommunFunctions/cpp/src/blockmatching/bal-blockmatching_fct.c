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
#include <math.h>
#include <string.h>
#include <time.h>

#ifndef WIN32
#include <unistd.h>
#include <sys/time.h>
#endif

#include <chunks.h>

#include <bal-behavior.h>
#include <bal-blockmatching-param.h>
#include <bal-blockmatching-param-tools.h>
#include <bal-blockmatching.h>
#include <bal-transformation-tools.h>
#include <bal-field-tools.h>


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

#define FILENAMELENGTH 128

static char *program = "blockMatching";

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
#ifndef WIN32  
  int print_time;
#endif
  /* misc
   */
  int print_parameters;


} local_parameter;

static void _PrintParam( FILE *f, local_parameter *par );


int blockmatching(
        char *floating_image,
        char *reference_image,
        char *result_image,
        char *initial_real_transformation,
        char *initial_voxel_transformation,
        char *initial_result_real_transformation,
        char *initial_result_voxel_transformation,
        char *result_real_transformation,
        char *result_voxel_transformation,
        int normalisation,
	bal_blockmatching_pyramidal_param param,
        int use_default_filename,
        char *command_line_file,
        char *log_file,
        int print_time,
        int print_parameters,
	int isDebug
        )
{
  char *auxptr;
  char auxstr[FILENAMELENGTH];
  local_parameter p;
  bal_image theFloatingImage;
  bal_image theReferenceImage;
  bal_image theResultImage;
  bal_transformation theInitTransformation;
  bal_transformation *initTransformation = (bal_transformation *)NULL;
  bal_transformation theTransformation;
  bal_transformation theResTransformation;
  bal_transformation *resTransformation = (bal_transformation *)NULL;

  p.floating_image = floating_image;
  p.reference_image = reference_image;
  p.result_image = result_image;
  p.initial_real_transformation = initial_real_transformation;
  p.initial_voxel_transformation = initial_voxel_transformation;
  p.initial_result_real_transformation = initial_result_real_transformation;
  p.initial_result_voxel_transformation = initial_result_voxel_transformation;
  p.result_real_transformation = result_real_transformation;
  p.result_voxel_transformation = result_voxel_transformation;
  p.normalisation = normalisation;
  p.param = param;
  p.use_default_filename = use_default_filename;
  p.command_line_file = command_line_file;
  p.log_file = log_file;

  #ifndef WIN32
  p.print_time = print_time;
  time_t t = time(NULL); 
  #endif
  
  p.print_parameters = print_parameters;


  /* printing parameters
   */
  if ( 0 && p.print_parameters ) {
    _PrintParam( stdout, &p );
    return( 0 );
  }



  /* writing stuff
     - write command line
     - open log file
   */
    if ( p.param.verbose > 0 ) {
    if ( p.log_file != NULL ) {
      if ( strlen( p.log_file ) == 4 && strcmp( p.log_file, "NULL" ) == 0 ) 
	p.param.verbosef = NULL;
      else if ( strlen( p.log_file ) == 6 && strcmp( p.log_file, "stderr" ) == 0 )
	p.param.verbosef = stderr;
      else if ( strlen( p.log_file ) == 6 && strcmp( p.log_file, "stdout" ) == 0 )
	p.param.verbosef = stdout;
      else {
	p.param.verbosef = fopen( p.log_file, "w" );
	if ( p.param.verbosef == NULL ) {
	  fprintf( stderr, "%s: unable to open '%s' for writing, switch to 'stderr'\n", 
		   program, p.log_file );
	  p.param.verbosef = stderr;
	}
      }
    }
    else if ( p.use_default_filename == 1 ) {
      sprintf( auxstr, "%s-%d-trace.txt", program, getpid() );
      p.param.verbosef = fopen( auxstr, "w" );
      if ( p.param.verbosef == NULL ) {
	fprintf( stderr, "%s: unable to open '%s' for writing, switch to 'stderr'\n", 
		 program, auxstr );
	p.param.verbosef = stderr;
      }
    }
#ifndef WIN32
    fprintf( p.param.verbosef, "%% %s\n", ctime( &t ) ); 
#endif
    BAL_SetVerboseFileInBalFieldTools( p.param.verbosef );
  }
  
  /* writing some stuff
   */
  if (p.param.verbosef != NULL) {
    BAL_PrintDefines( p.param.verbosef );
  }

  if ( p.print_parameters ) {
    if (p.param.verbosef != NULL) _PrintParam( p.param.verbosef, &p );
    else _PrintParam( stdout, &p );
  }


  


   /***************************************************
   *
   * 
   *
   ***************************************************/
 
  /* reading reference image
   */
  if ( p.param.verbosef != NULL ) {
    fprintf( p.param.verbosef, "\tReading reference image '%s'\n", p.reference_image );
  }
  if ( BAL_ReadImage( &theReferenceImage, p.reference_image, p.normalisation ) != 1 ) {
    if ( p.param.verbosef != NULL && p.param.verbosef != stderr && p.param.verbosef != stdout )
      fclose( p.param.verbosef );
    fprintf( stderr, "%s: can not read '%s'\n", program, p.reference_image );
    exit( -1 );
  }

  /* reading floating image
   */
  if ( p.param.verbosef != NULL ) {
    fprintf( p.param.verbosef, "\tReading floating image '%s'\n", p.floating_image );
  }
  if ( BAL_ReadImage( &theFloatingImage, p.floating_image, p.normalisation ) != 1 ) {    
    BAL_FreeImage( &theReferenceImage );
    if ( p.param.verbosef != NULL && p.param.verbosef != stderr && p.param.verbosef != stdout )
      fclose( p.param.verbosef );
    fprintf( stderr, "%s: can not read '%s'\n", program, p.floating_image );
    exit( -1 );
  }
  
  /* adjust parameters
     image may be 2D
  */
  BAL_AdjustBlockMatchingPyramidalParameters( &theReferenceImage, &theFloatingImage, &(p.param) );





  /***************************************************
   *
   * 
   *
   ***************************************************/

  /* allocating result transformation
     it is initialized to the identity
   */
  /* reading  or compute initial result transformation
   */
  if ( p.initial_result_real_transformation!= NULL ) {
    
    if ( p.param.verbosef != NULL ) {
      fprintf( p.param.verbosef, "\tReading initial 'real' result transformation '%s'\n", 
	       p.initial_result_real_transformation );
    }

    if ( BAL_ReadTransformation( &theTransformation, p.initial_result_real_transformation ) != 1 ) {
      BAL_FreeImage( &theFloatingImage );
      BAL_FreeImage( &theReferenceImage );
      fprintf( stderr, "%s: unable to read 'real' result transformation '%s'\n", 
	       program, p.initial_result_real_transformation );
      exit( -1 );
    }
    theTransformation.transformation_unit = REAL_UNIT;
    
  }
  else if ( p.initial_result_voxel_transformation != NULL ) {
    
    if ( p.param.verbosef != NULL ) {
      fprintf( p.param.verbosef, "\tReading initial 'voxel' result transformation '%s'\n", 
	       p.initial_result_voxel_transformation );
    }

    if ( BAL_ReadTransformation( &theTransformation, p.initial_result_voxel_transformation ) != 1 ) {
      BAL_FreeImage( &theFloatingImage );
      BAL_FreeImage( &theReferenceImage );
      fprintf( stderr, "%s: unable to read 'voxel' result transformation '%s'\n", 
	       program, p.initial_result_voxel_transformation );
      exit( -1 );
    }
    theTransformation.transformation_unit = VOXEL_UNIT;
    
    if ( BAL_ChangeTransformationToRealUnit( &theReferenceImage, &theFloatingImage, 
					     &theTransformation, &theTransformation ) != 1 ) {
      BAL_FreeTransformation( &theTransformation );
      BAL_FreeImage( &theFloatingImage );
      BAL_FreeImage( &theReferenceImage );
      fprintf( stderr, "%s: unable to convert 'voxel' transformation '%s' into the 'real' world\n", 
	       program, p.initial_result_voxel_transformation );
      exit( -1 );
    }

  }
  else {

    if ( BAL_AllocTransformation( &theTransformation, p.param.transformation_type, &theReferenceImage ) != 1 ) {
      BAL_FreeImage( &theFloatingImage );
      BAL_FreeImage( &theReferenceImage );
      if ( p.param.verbosef != NULL && p.param.verbosef != stderr && p.param.verbosef != stdout )
	fclose( p.param.verbosef );
      fprintf( stderr, "%s: can not allocate initial result transformation\n", program );
      exit( -1 );
    }
    
  }













  /* reading  or compute initial transformation
   */
  if ( p.initial_real_transformation != NULL ) {
    
    if ( p.param.verbosef != NULL ) {
      fprintf( p.param.verbosef, "\tReading initial 'real' transformation '%s'\n", p.initial_real_transformation );
    }

    if ( BAL_ReadTransformation( &theInitTransformation, p.initial_real_transformation ) != 1 ) {
      BAL_FreeTransformation( &theTransformation );
      BAL_FreeImage( &theFloatingImage );
      BAL_FreeImage( &theReferenceImage );
      fprintf( stderr, "%s: unable to read 'real' transformation '%s'\n", program, p.initial_real_transformation );
      exit( -1 );
    }
    theInitTransformation.transformation_unit = REAL_UNIT;
    initTransformation = &theInitTransformation;

  }
  else if ( p.initial_voxel_transformation != NULL ) {
    
    if ( p.param.verbosef != NULL ) {
      fprintf( p.param.verbosef, "\tReading initial 'voxel' transformation '%s'\n", p.initial_voxel_transformation );
    }

    if ( BAL_ReadTransformation( &theInitTransformation, p.initial_voxel_transformation ) != 1 ) {
      BAL_FreeTransformation( &theTransformation );
      BAL_FreeImage( &theFloatingImage );
      BAL_FreeImage( &theReferenceImage );
      fprintf( stderr, "%s: unable to read 'voxel' transformation '%s'\n", program, p.initial_voxel_transformation );
      exit( -1 );
    }
    theInitTransformation.transformation_unit = VOXEL_UNIT;
    
    if ( isDebug )
      BAL_PrintTransformation( stderr, &theInitTransformation, "p.initial_voxel_transformation" );

    if ( BAL_ChangeTransformationToRealUnit( &theReferenceImage, &theFloatingImage, 
					     &theInitTransformation, &theInitTransformation ) != 1 ) {
      BAL_FreeTransformation( &theInitTransformation );
      BAL_FreeTransformation( &theTransformation );
      BAL_FreeImage( &theFloatingImage );
      BAL_FreeImage( &theReferenceImage );
      fprintf( stderr, "%s: unable to convert 'voxel' transformation '%s' into the 'real' world\n", 
	       program, p.initial_voxel_transformation );
      exit( -1 );
    }

    initTransformation = &theInitTransformation;

  }
  else if ( theReferenceImage.ncols != theFloatingImage.ncols
	    || theReferenceImage.nrows != theFloatingImage.nrows
	    || theReferenceImage.nplanes != theFloatingImage.nplanes
	    || fabs( theReferenceImage.vx - theFloatingImage.vx )/theReferenceImage.vx > 0.001
	    || fabs( theReferenceImage.vy - theFloatingImage.vy )/theReferenceImage.vy > 0.001
	    || fabs( theReferenceImage.vz - theFloatingImage.vz )/theReferenceImage.vz > 0.001
	    ) {

    if ( p.param.verbosef != NULL ) {
      fprintf( p.param.verbosef, "\tComputing '%s' - '%s' transformation\n", 
	       theReferenceImage.name, theFloatingImage.name );
    }  

    if ( BAL_AllocTransformation( &theInitTransformation, AFFINE_3D, (bal_image *)NULL ) != 1 ) {
      BAL_FreeTransformation( &theTransformation );
      BAL_FreeImage( &theFloatingImage );
      BAL_FreeImage( &theReferenceImage );
      if ( p.param.verbosef != NULL && p.param.verbosef != stderr && p.param.verbosef != stdout )
	fclose( p.param.verbosef );
      fprintf( stderr, "%s: can not allocate initial transformation\n", program );
      exit( -1 );
    }
    
    if ( BAL_ComputeImageToImageTransformation( &theReferenceImage, &theFloatingImage, &theInitTransformation ) != 1 ) {
      BAL_FreeTransformation( &theInitTransformation );
      BAL_FreeTransformation( &theTransformation );
      BAL_FreeImage( &theFloatingImage );
      BAL_FreeImage( &theReferenceImage );
      if ( p.param.verbosef != NULL && p.param.verbosef != stderr && p.param.verbosef != stdout )
	fclose( p.param.verbosef );
      fprintf( stderr, "%s: can not compute initial transformation\n", program );
      exit( -1 );
    }

    initTransformation = &theInitTransformation;
    
  }




  /***************************************************
   *
   * matching
   *
   ***************************************************/
  
  /* do nothing if there is no iterations per level
   */

  if ( p.param.max_iterations.lowest > 0 || p.param.max_iterations.highest > 0 ) {
    
    if ( isDebug ) {
      fprintf( stderr, "========== before matching ==========\n" );
      BAL_PrintTransformation( stderr,  initTransformation, "initial transformation" );
      BAL_PrintTransformation( stderr,  &theTransformation, "transformation to be computed" );
      fprintf( stderr, "=====================================\n" );    }
    
    if ( BAL_PyramidalBlockMatching( &theReferenceImage, &theFloatingImage,
				     &theTransformation, initTransformation,
				     &(p.param) ) != 1 ) {
      BAL_FreeTransformation( &theInitTransformation );
      BAL_FreeTransformation( &theTransformation );
      BAL_FreeImage( &theFloatingImage );
      BAL_FreeImage( &theReferenceImage );
      if ( p.param.verbosef != NULL && p.param.verbosef != stderr && p.param.verbosef != stdout )
	fclose( p.param.verbosef );
      fprintf( stderr, "%s: unable to register the images \n", program );
      exit( -1 );
    }
    
    if ( isDebug ) {
      fprintf( stderr, "========== after matching ===========\n" );
      BAL_PrintTransformation( stderr,  &theTransformation, "computed transformation" );
      fprintf( stderr, "=====================================\n" );
    }

  }






  /***************************************************
   *
   * resample floating image
   *
   ***************************************************/
  
  /* images are to be re-read again
     - re-read the floating image
  */
  if ( p.normalisation ) {
    BAL_FreeImage( &theFloatingImage );
    if ( BAL_ReadImage( &theFloatingImage, p.floating_image, 0 ) != 1 ) {    
      BAL_FreeTransformation( &theInitTransformation );
      BAL_FreeTransformation( &theTransformation );
      BAL_FreeImage( &theResultImage );
      BAL_FreeImage( &theReferenceImage );
      if ( p.param.verbosef != NULL && p.param.verbosef != stderr && p.param.verbosef != stdout )
	fclose( p.param.verbosef );
      fprintf( stderr, "%s: can not read '%s' (resampling step)\n", program, p.floating_image );
      exit( -1 );
    }
  }



  /* allocate a result image with the same type than the floating one
     and with the geometry of the reference one
  */
  sprintf( auxstr, "%s-%d-result-image.hdr", program, getpid() );
  if ( BAL_InitAllocImage( &theResultImage, auxstr, 
			   theReferenceImage.ncols, theReferenceImage.nrows, 
			   theReferenceImage.nplanes, theReferenceImage.vdim, theFloatingImage.type ) == -1 ) {
    BAL_FreeTransformation( &theInitTransformation );
    BAL_FreeTransformation( &theTransformation );
    BAL_FreeImage( &theFloatingImage );
    BAL_FreeImage( &theReferenceImage );
    if ( p.param.verbosef != NULL && p.param.verbosef != stderr && p.param.verbosef != stdout )
      fclose( p.param.verbosef );
    fprintf( stderr, "%s: can not allocate result image (resampling step)\n", program );
    exit( -1 );
  }
  theResultImage.vx = theReferenceImage.vx;
  theResultImage.vy = theReferenceImage.vy;
  theResultImage.vz = theReferenceImage.vz;


  
  /* compose init o restrsf
   */
  if ( initTransformation != (bal_transformation *)NULL ) {

    if ( BAL_AllocTransformationComposition( &theResTransformation,
					     initTransformation, 
					     &theTransformation,
					     &theReferenceImage ) != 1 ) {
      BAL_FreeImage( &theResultImage );
      BAL_FreeTransformation( &theInitTransformation );
      BAL_FreeTransformation( &theTransformation );
      BAL_FreeImage( &theFloatingImage );
      BAL_FreeImage( &theReferenceImage );
      if ( p.param.verbosef != NULL && p.param.verbosef != stderr && p.param.verbosef != stdout )
	fclose( p.param.verbosef );
      fprintf( stderr, "%s: can not allocate composition transformation (resampling step)\n", program );
      exit( -1 );
    }

    if ( BAL_TransformationComposition( &theResTransformation,
					initTransformation, 
					&theTransformation  ) != 1 ) {
      BAL_FreeTransformation( &theResTransformation );
      BAL_FreeImage( &theResultImage );
      BAL_FreeTransformation( &theInitTransformation );
      BAL_FreeTransformation( &theTransformation );
      BAL_FreeImage( &theFloatingImage );
      BAL_FreeImage( &theReferenceImage );
      if ( p.param.verbosef != NULL && p.param.verbosef != stderr && p.param.verbosef != stdout )
	fclose( p.param.verbosef );
      fprintf( stderr, "%s: can not compute composition transformation (resampling step)\n", program );
      exit( -1 );
    }
    resTransformation = &theResTransformation;

  }
  else {
    resTransformation = &theTransformation;
  }



  /* the initial transformation is no longer required
   */  
  if ( initTransformation != (bal_transformation *)NULL ) {
    BAL_FreeTransformation( &theInitTransformation );
    initTransformation = (bal_transformation *)NULL;
  }

  
  
  /* resample floating image
   */
  if ( BAL_ResampleImage( &theFloatingImage, &theResultImage, resTransformation, LINEAR ) != 1 ) {
    BAL_FreeTransformation( &theResTransformation );
    BAL_FreeImage( &theResultImage );
    BAL_FreeTransformation( &theTransformation );
    BAL_FreeImage( &theFloatingImage );
    BAL_FreeImage( &theReferenceImage );
    if ( p.param.verbosef != NULL && p.param.verbosef != stderr && p.param.verbosef != stdout )
      fclose( p.param.verbosef );
    fprintf( stderr, "%s: can not resample floating image (resampling step)\n", program );
    exit( -1 );
  }


 
  /***************************************************
   *
   * writing the (resampling) transformation
   * it is the composition "initial o computed"
   *
   ***************************************************/

  if ( isDebug )
    BAL_PrintTransformation( stderr,  resTransformation, "resampling transformation" );

  /* writing the "real" transformation
   */

  auxptr = (char*)NULL;

  if ( p.result_real_transformation != NULL ) {
    auxptr = p.result_real_transformation;
  } 
  else if ( p.use_default_filename == 1 ) {
    sprintf( auxstr, "%s-%d-result-transformation.trsf", program, getpid() );
    auxptr = auxstr;
  }

  if ( auxptr != (char*)NULL ) {
    if ( BAL_WriteTransformation( resTransformation, auxptr ) != 1 ) {
      BAL_FreeTransformation( &theResTransformation );
      BAL_FreeImage( &theResultImage );
      BAL_FreeTransformation( &theTransformation );
      BAL_FreeImage( &theFloatingImage );
      BAL_FreeImage( &theReferenceImage );
      if ( p.param.verbosef != NULL && p.param.verbosef != stderr && p.param.verbosef != stdout )
	fclose( p.param.verbosef );
      fprintf( stderr, "%s: unable to write real result transformation '%s'\n", program, auxptr );
      exit( -1 );
    }
  }
 
  /* writing the "voxel" transformation
   */

  if ( p.result_voxel_transformation != NULL ) {
    if ( BAL_ChangeTransformationToVoxelUnit( &theFloatingImage, &theReferenceImage, 
					      resTransformation, resTransformation ) != 1 ) {
      BAL_FreeTransformation( &theResTransformation );
      BAL_FreeImage( &theResultImage );
      BAL_FreeTransformation( &theTransformation );
      BAL_FreeImage( &theFloatingImage );
      BAL_FreeImage( &theReferenceImage );
      fprintf( stderr, "%s: unable to convert 'real' transformation into the 'voxel' world\n", 
	       program );
      exit( -1 );
    }

    auxptr = p.result_voxel_transformation;

    if ( BAL_WriteTransformation( resTransformation, auxptr ) != 1 ) {
      BAL_FreeTransformation( &theResTransformation );
      BAL_FreeImage( &theResultImage );
      BAL_FreeTransformation( &theTransformation );
      BAL_FreeImage( &theFloatingImage );
      BAL_FreeImage( &theReferenceImage );
      if ( p.param.verbosef != NULL && p.param.verbosef != stderr && p.param.verbosef != stdout )
	fclose( p.param.verbosef );
      fprintf( stderr, "%s: unable to write voxel result transformation '%s'\n", program, auxptr );
      exit( -1 );
    }
  }


  

  /***************************************************
   *
   * freeing stuff
   *
   ***************************************************/
 
  /* only allocated when an initial transformation is given
   */
  if ( initTransformation != (bal_transformation *)NULL ) {
    BAL_FreeTransformation( &theResTransformation );
  }
  BAL_FreeTransformation( &theTransformation );
  BAL_FreeImage( &theFloatingImage );
  BAL_FreeImage( &theReferenceImage );



  /***************************************************
   *
   * writing result image
   *
   ***************************************************/

  if ( p.result_image != NULL ) {
    auxptr = p.result_image;
    } 
  else {
    if ( p.use_default_filename == 1 ) {
      sprintf( auxstr, "%s-%d-result-image.hdr", program, getpid() );
      auxptr = auxstr;
    }
    else {
      BAL_FreeImage( &theResultImage );
      if ( p.param.verbosef != NULL && p.param.verbosef != stderr && p.param.verbosef != stdout )
	fclose( p.param.verbosef );
      fprintf( stderr, "%s: no name for result image (resampling step)\n", program );
      exit( -1 );
    }
  }
  
  if ( BAL_WriteImage( &theResultImage, auxptr ) != 1 ) {
    BAL_FreeImage( &theResultImage );
    if ( p.param.verbosef != NULL && p.param.verbosef != stderr && p.param.verbosef != stdout )
      fclose( p.param.verbosef );
    fprintf( stderr, "%s: unable to write result image (resampling step)\n", program );
    exit( -1 );
  }



  /* free result image
   */
  BAL_FreeImage( &theResultImage );
  


  /* writing stuff
     - close log file
   */
  if ( p.param.verbosef != NULL && p.param.verbosef != stderr && p.param.verbosef != stdout )
    fclose( p.param.verbosef );



  return 0;
}


static void _PrintParam( FILE *f, local_parameter *p ) 
{
  fprintf( f, "--- file names\n" );
  
  fprintf( f, "p->floating_image = " );
  if ( p->floating_image != NULL ) fprintf( f, "'%s'\n", p->floating_image );
  else fprintf( f, "NULL\n" );
  fprintf( f, "p->reference_image = " );
  if ( p->reference_image != NULL ) fprintf( f, "'%s'\n", p->reference_image );
  else fprintf( f, "NULL\n" );
  fprintf( f, "p->result_image = " );
  if ( p->result_image != NULL ) fprintf( f, "'%s'\n", p->result_image );
  else fprintf( f, "NULL\n" );
  
  fprintf( f, "p->initial_real_transformation = " );
  if ( p->initial_real_transformation != NULL ) fprintf( f, "'%s'\n", p->initial_real_transformation );
  else fprintf( f, "NULL\n" );
  fprintf( f, "p->initial_voxel_transformation = " );
  if ( p->initial_real_transformation != NULL ) fprintf( f, "'%s'\n", p->initial_voxel_transformation );
  else fprintf( f, "NULL\n" );
  fprintf( f, "p->result_real_transformation = " );
  if ( p->result_real_transformation != NULL ) fprintf( f, "'%s'\n", p->result_real_transformation );
  else fprintf( f, "NULL\n" );
  fprintf( f, "p->result_voxel_transformation = " );
  if ( p->result_voxel_transformation != NULL ) fprintf( f, "'%s'\n", p->result_voxel_transformation );
  else fprintf( f, "NULL\n" );

  fprintf( f, "--- pre-processing\n" );
  fprintf( f, "p->normalisation = %d\n", p->normalisation );

  fprintf( f, "--- writing stuff\n" );
  fprintf( f, "p->use_default_filename = %d\n", p->use_default_filename );
  fprintf( f, "p->command_line_file = " );
  if ( p->command_line_file != NULL ) fprintf( f, "'%s'\n", p->command_line_file );
  else fprintf( f, "NULL\n" );
  fprintf( f, "p->log_file = " );
  if ( p->log_file != NULL ) fprintf( f, "'%s'\n", p->log_file );
  else fprintf( f, "NULL\n" );
  
  fprintf( f, "--- misc\n" );
  fprintf( f, "p->print_parameters =  %d\n", p->print_parameters );
#ifndef WIN32
  fprintf( f, "p->print_time =  %d\n", p->print_time );
#endif
  fprintf( f, "\n" );
  fprintf( f, "--- parameters for hierarchical block matching\n" );
  BAL_PrintBlockMatchingPyramidalParameters( f, &(p->param) );
}
