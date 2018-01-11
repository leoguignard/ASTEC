/*************************************************************************
 * changeMultipleTrsfs.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mar 21 jan 2014 22:43:26 CET
 *
 * ADDITIONS, CHANGES
 *
 */


#include <sys/time.h> /* gettimeofday() */
#include <time.h> /* clock() */
#include <string.h>

#include <string-tools.h>

#include <bal-transformation-list-tools.h>
#include <bal-transformation-averaging.h>

static int _verbose_ = 1;
static int _debug_ = 0;





typedef enum {
  _MEMORY_,
  _STREAMING_
} typeComputation;


typedef struct local_parameter {

  char *thetransformation_format;
  int firstindex;
  int lastindex;
  int refindex;

  char *theimage_name;

  char *restransformation_format;

  char *resimage_name;

  int isotropic;
  int margin;

  enumTypeTransfo transformation_type;

  typeComputation computation;
  int print_time;

} local_parameter;





/*------- Definition des fonctions statiques ----------*/
static void _ErrorParse( char *str, int flag );
static void _Parse( int argc, char *argv[], local_parameter *p );
static void _InitParam( local_parameter *par );
static double _GetTime();
static double _GetClock();




static char *program = NULL;

static char *usage = "[format-in] -f[irst] %d -l[ast] %d [-r[ef|eference] %d]\n\
 [-template %s]\n\
 [[-res] format-out] [-res-template %s] [-iso] [-margin %d]\n\
 [-transformation-type|-transformation|-trsf-type %s]\n\
 [-streaming | -memory]\n\
 [-v] [-help]";

static char *detail = "\
\n\
given a list of transformation towards a reference (of the same size than the\n\
template), compute a new template that contains all transformed images as well\n\
as the new transformations\n\
\n\
[format-in]    # format 'a la printf' of transformations to be processed\n\
               # must contain one '%d'\n\
               # depicts transformations of the form T_{i<-ref}\n\
[-first %d]    # first value of the index in the format\n\
[-last %d]     # last value of the index in the format\n\
[-reference %d] # index of the reference transformation\n\
                # will only be translated\n\
[-template %s] # template image corresponding to the input images\n\
[-res] format-out  # format 'a la printf' for output transformations\n\
               # will allow to resample input image into the resulting\n\
               # template (thus still of the form T_{i<-ref})\n\
[-res-template %s] # template image corresponding to the output images\n\
[-iso]         # make voxels isotropic\n\
[-margin %d]   # add a margin (in voxels)\n\
[-transformation-type|-transformation|-trsf-type %s] # transformation type\n\
  translation2D, translation3D, translation-scaling2D, translation-scaling3D,\n\
  rigid2D, rigid3D, rigid, similitude2D, similitude3D, similitude,\n\
  affine2D, affine3D, affine, vectorfield2D, vectorfield3D, vectorfield, vector\n\
[-streaming] # computation is done by reading one transformation after the other\n\
[-memory]    # computation is done by loading all transformation in memory\n\
 -v : mode verbose\n\
\n";







int main( int argc, char *argv[] )
{
  local_parameter p;

  stringList theTransformationFileList;
  bal_transformationList theTransformationList; 
  bal_image theTemplate;

  int ntransformations = 0;

  stringList resTransformationFileList;
  bal_transformationList resTransformationList; 
  bal_image resTemplate;

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
  

  
  ntransformations = p.lastindex - p.firstindex + 1;





  /***************************************************
   *
   * inputs
   *
   ***************************************************/


  /* transformation file names
     string[j][i] is printf( format, i, j )
     since transformations are denoted floating_on_reference and allows to resample
     the floating image onto the reference
     string[j][i] refers to the transformation T_{i<-j}
   */
  initStringList( &theTransformationFileList );

  if ( buildStringListFromFormat( p.thetransformation_format, 
				  p.firstindex, p.lastindex, 
				  &theTransformationFileList ) != 1 ) {
    _ErrorParse( "unable to build input transformation list from format\n", 0);
  }
  
  if ( _debug_ >= 3 ) 
    printStringList( stderr, &theTransformationFileList, "Input transformations" );





  /* reading transformations
   */
  BAL_InitTransformationList( &theTransformationList );

  if ( BAL_AllocTransformationList( &theTransformationList, 
				     ntransformations ) != 1 ) {
    freeStringList( &theTransformationFileList );
    _ErrorParse( "error when allocating input transformations\n", 0);
  }

  
  if ( BAL_ReadTransformationList( &theTransformationList, &theTransformationFileList ) != 1 ) {
    BAL_FreeTransformationList( &theTransformationList );
    freeStringList( &theTransformationFileList );
    _ErrorParse( "error when reading input transformations\n", 0);
  }

  freeStringList( &theTransformationFileList );

  

  /* reading input template
   */
  BAL_InitImage( &theTemplate, NULL, 0, 0, 0, 0, UCHAR );
  if ( BAL_ReadImage( &theTemplate, p.theimage_name, 0 ) != 1 ) {
    BAL_FreeTransformationList( &theTransformationList );
    _ErrorParse( "error when reading input template\n", 0);
  }
  
  

 

  /***************************************************
   *
   * outputs
   *
   ***************************************************/
  

  /* initializing result transformations
     should we check whether read transformations are compatible with
     given type ?
  */
  BAL_InitTransformationList( &resTransformationList );

  switch ( p.transformation_type ) {
  default :
    BAL_FreeImage( &theTemplate );
    BAL_FreeTransformationList( &theTransformationList );
    _ErrorParse( "such transformation type not handled yet (or no type was given)\n", 0);
    
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
    
    if ( BAL_FullAllocTransformationList( &resTransformationList, 
					  ntransformations,
					  p.transformation_type, (bal_image *)NULL ) != 1 ) {
      BAL_FreeImage( &theTemplate );
      BAL_FreeTransformationList( &theTransformationList );
      _ErrorParse( "unable to allocate result transformations (matrices)\n", 0);
    }
    break;
    
  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
    
    /* should provide a template
     */
    BAL_FreeImage( &theTemplate );
    BAL_FreeTransformationList( &theTransformationList );
    _ErrorParse( "such transformation type not handled yet (vector field)\n", 0);
  }

  /* result template
   */
  BAL_InitImage( &resTemplate, NULL, 0, 0, 0, 0, UCHAR );





  /***************************************************
   *
   * 
   *
   ***************************************************/
  if ( p.refindex < p.firstindex || p.lastindex < p.refindex )
    p.refindex = p.firstindex;

  if ( BAL_ChangeTransformationList( &theTransformationList, &theTemplate,
				     &resTransformationList, &resTemplate,
				     p.refindex - p.firstindex,
				     p.margin, p.isotropic ) != 0 ) {
    
    BAL_FreeTransformationList( &resTransformationList );
    BAL_FreeImage( &theTemplate );
    BAL_FreeTransformationList( &theTransformationList );
    _ErrorParse( "unable to recompute transformations\n", 0);
  }

  BAL_FreeImage( &theTemplate );
  BAL_FreeTransformationList( &theTransformationList );




  
   /* result transformation file names
   */

  initStringList( &resTransformationFileList );

  if ( p.restransformation_format != (char*)NULL ) {
    if ( buildStringListFromFormat( p.restransformation_format, p.firstindex, p.lastindex, 
				    &resTransformationFileList ) != 1 ) {
      BAL_FreeTransformationList( &resTransformationList );
      _ErrorParse( "unable to build output transformation list from format\n", 0);
    }
  }
  else {
    BAL_FreeTransformationList( &resTransformationList );
    _ErrorParse( "no format for the result transformations\n", 0);
  }

  
  if ( _debug_ >= 3 ) 
    printStringList( stderr, &resTransformationFileList, "Result transformations" );

  

  /* writing results
   */


  if ( BAL_WriteTransformationList( &resTransformationList, &resTransformationFileList ) != 1 ) {
    freeStringList( &resTransformationFileList );
    BAL_FreeTransformationList( &resTransformationList );
    _ErrorParse( "unable to write result transformations\n", 0);
  }



  freeStringList( &resTransformationFileList );
  BAL_FreeTransformationList( &resTransformationList );


  /* ... */
  if ( p.resimage_name != (char*)NULL ) {
    if ( BAL_AllocImage( &resTemplate ) != 1 ) {
      _ErrorParse( "unable to allocate result template image\n", 0);
    }
    if ( BAL_WriteImage( &resTemplate, p.resimage_name ) != 1 ) {
      _ErrorParse( "unable to write result template image\n", 0);
    }
    BAL_FreeImage( &resTemplate );
  }







  

  
  time_exit = _GetTime();
  clock_exit = _GetClock();

  if (  p.print_time ) {
    fprintf( stderr, "%s: elapsed time = %f\n", program, time_exit - time_init );
    fprintf( stderr, "%s: elapsed time = %f\n", program, clock_exit - clock_init );
  }

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

      /* file names
       */

      if ( strcmp ( argv[i], "-format" ) == 0 ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -format...\n", 0 );
	if ( p->thetransformation_format != (char*)NULL || inputisread > 0 ) 
	  _ErrorParse( "parsing -format: input has already been parsed ...\n", 0 );
	p->thetransformation_format = argv[i];
	inputisread = 1;
      }

      else if ( (strcmp ( argv[i], "-f" ) == 0 && argv[i][2] == '\0') 
		|| (strcmp ( argv[i], "-first" ) == 0 && argv[i][6] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -first ...\n", 0 );
	status = sscanf( argv[i], "%d", &(p->firstindex) );
	if ( status <= 0 ) _ErrorParse( "parsing -first ...", 0 );
      }
      else if ( (strcmp ( argv[i], "-l" ) == 0 && argv[i][2] == '\0') 
		|| (strcmp ( argv[i], "-last" ) == 0 && argv[i][5] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -last ...\n", 0 );
	status = sscanf( argv[i], "%d", &(p->lastindex) );
	if ( status <= 0 ) _ErrorParse( "parsing -last ...", 0 );
      }
      else if ( (strcmp ( argv[i], "-r" ) == 0 && argv[i][2] == '\0') 
		|| (strcmp ( argv[i], "-ref" ) == 0 && argv[i][4] == '\0') 
		|| (strcmp ( argv[i], "-reference" ) == 0) ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -reference ...\n", 0 );
	status = sscanf( argv[i], "%d", &(p->refindex) );
	if ( status <= 0 ) _ErrorParse( "parsing -reference ...", 0 );
      }
      
      else if ( strcmp ( argv[i], "-template") == 0
		|| (strcmp ( argv[i], "-t") == 0 && argv[i][2] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -template", 0 );
	p->theimage_name = argv[i];
      }

      else if ( strcmp ( argv[i], "-res" ) == 0  && argv[i][4] == '\0' ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -res...\n", 0 );
	if ( p->restransformation_format != (char*)NULL || outputisread > 0 ) 
	  _ErrorParse( "parsing -res: output has already been parsed ...\n", 0 );
	p->restransformation_format = argv[i];
	outputisread = 1;
      }

      else if ( strcmp ( argv[i], "-result-template") == 0
		|| strcmp ( argv[i], "-res-template") == 0
		|| (strcmp ( argv[i], "-res-t") == 0 && argv[i][6] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -res-template", 0 );
	p->resimage_name = argv[i];
      }


      /* transformation type
       */
      else if ( strcmp ( argv[i], "-transformation-type" ) == 0 
		|| strcmp ( argv[i], "-transformation" ) == 0
		|| strcmp ( argv[i], "-trsf-type" ) == 0 ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "-transformation-type", 0 );
	if ( strcmp ( argv[i], "translation2D" ) == 0 ) {
	  p->transformation_type = TRANSLATION_2D;
	}
	else if ( strcmp ( argv[i], "translation3D" ) == 0 ) {
	  p->transformation_type = TRANSLATION_3D;
	}
	else if ( strcmp ( argv[i], "translation" ) == 0 && argv[i][11] == '\0') {
	  p->transformation_type = TRANSLATION_3D;
	}
	else if ( strcmp ( argv[i], "translation-scaling2D" ) == 0 ) {
	  p->transformation_type = TRANSLATION_SCALING_2D;
	}
	else if ( strcmp ( argv[i], "translation-scaling3D" ) == 0 ) {
	  p->transformation_type = TRANSLATION_SCALING_3D;
	}
	else if ( strcmp ( argv[i], "rigid2D" ) == 0 ) {
	  p->transformation_type = RIGID_2D;
	}
	else if ( strcmp ( argv[i], "rigid3D" ) == 0 ) {
	  p->transformation_type = RIGID_3D;
	}
	else if ( (strcmp ( argv[i], "rigid" ) == 0 && argv[i][5] == '\0') ) {
	  p->transformation_type = RIGID_3D;
	}
	else if ( strcmp ( argv[i], "similitude2D" ) == 0 ) {
	  p->transformation_type = SIMILITUDE_2D;
	}
	else if ( strcmp ( argv[i], "similitude3D" ) == 0 ) {
	  p->transformation_type = SIMILITUDE_3D;
	}
	else if ( strcmp ( argv[i], "similitude" ) == 0 ) {
	  p->transformation_type = SIMILITUDE_3D;
	}
	else if ( strcmp ( argv[i], "affine2D" ) == 0 ) {
	  p->transformation_type = AFFINE_2D;
	}
	else if ( strcmp ( argv[i], "affine3D" ) == 0 ) {
	  p->transformation_type = AFFINE_3D;
	}
	else if ( strcmp ( argv[i], "affine" ) == 0 ) {
	  p->transformation_type = AFFINE_3D;
	}
	/*
	  else if ( strcmp ( argv[i], "spline" ) == 0 ) {
	  p->transformation_type = SPLINE;
	  }
	*/
	else if ( strcmp ( argv[i], "vectorfield" ) == 0 
		  || (strcmp ( argv[i], "vector" ) == 0 && argv[i][6] == '\0') ) {
	  p->transformation_type = VECTORFIELD_3D;
	}
	else if ( strcmp ( argv[i], "vectorfield3D" ) == 0 
		  || strcmp ( argv[i], "vector3D" ) == 0 ) {
	  p->transformation_type = VECTORFIELD_3D;
	}
	else if ( strcmp ( argv[i], "vectorfield2D" ) == 0 
		  || strcmp ( argv[i], "vector2D" ) == 0 ) {
	  p->transformation_type = VECTORFIELD_2D;
	}
	else {
	  fprintf( stderr, "unknown transformation type: '%s'\n", argv[i] );
	  _ErrorParse( "-transformation-type", 0 );
	}
      }

      

      /* computation tuning
       */
      else if ( (strcmp ( argv[i], "-iso" ) == 0 && argv[i][4] == '\0') 
		|| strcmp ( argv[i], "-isotropic" ) == 0 ) {
	p->isotropic = 1;
      }
      
      else if ( strcmp ( argv[i], "-margin" ) == 0 ) {
	if ( i >= argc) _ErrorParse( "parsing -margin ...\n", 0 );
	status = sscanf( argv[i], "%d", &(p->margin) );
	if ( status <= 0 ) _ErrorParse( "parsing -margin ...", 0 );
      }

      

      /* kind of computation
       */
      else if ( (strcmp ( argv[i], "-memory" ) == 0 && argv[i][7] == '\0') ) {
	p->computation = _MEMORY_;
      }
      else if ( (strcmp ( argv[i], "-streaming" ) == 0 && argv[i][10] == '\0') ) {
	p->computation = _STREAMING_;
      }

      /* general arguments*/
       else if ( (strcmp ( argv[i], "-time" ) == 0 && argv[i][5] == '\0') ) {
	p->print_time = 1;
      }
      else if ( (strcmp ( argv[i], "-notime" ) == 0 && argv[i][7] == '\0')  
		|| (strcmp ( argv[i], "-no-time" ) == 0 && argv[i][8] == '\0') ) {
	p->print_time = 0;
      }

      else if ( (strcmp ( argv[i], "-verbose") == 0 && argv[i][8] == '\0')
		|| (strcmp ( argv[i], "-v") == 0 && argv[i][2] == '\0') ) {
	if ( _verbose_ <= 0 ) _verbose_ = 1;
	else _verbose_ ++;
	BAL_IncrementVerboseInBalTransformationAveraging(  );
	BAL_IncrementVerboseInBalTransformation(  );
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
	p->thetransformation_format = argv[i];
	inputisread = 1;
      }
      else if ( outputisread == 0 ) {
	p->restransformation_format = argv[i];
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
  p->thetransformation_format = (char*)NULL;
  p->firstindex = 0;
  p->lastindex = 0;
  p->refindex = 0;

  p->theimage_name = (char*)NULL;

  p->restransformation_format = (char*)NULL;

  p->transformation_type = UNDEF_TRANSFORMATION;

  p->resimage_name = (char*)NULL;

  p->isotropic = 0;

  p->margin = 0;

  p->computation = _MEMORY_;

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

