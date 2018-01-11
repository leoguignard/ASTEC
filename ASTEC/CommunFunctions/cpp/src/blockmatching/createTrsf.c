/*************************************************************************
 * createRandomTrsf.c -
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

#include <bal-transformation-tools.h>
#include <bal-createRandomTrsf.h>

typedef enum {
  _RANDOM_,
  _IDENTITY_
} enumValueTransfo;


typedef struct local_parameter {
  char *restrsf_name;
  char *template_name;

  bal_doublePoint fixedpoint;

  enumTypeTransfo transformation_type;
  enumValueTransfo transformation_value;
  int print;
} local_parameter;


static char *program = NULL;

static char *usage = "%s\n\
 [-transformation-type|-transformation|-trsf-type %s] [-print | -no-print]\n\
 [-random | -identity]\n\
 [-template %s] [-fixedpoint %lf %lf [%lf]]\n";

static char *detail = "\
[-transformation-type|-transformation|-trsf-type %s] # transformation type\n\
  translation2D, translation3D, translation-scaling2D, translation-scaling3D,\n\
  rigid2D, rigid3D, rigid, similitude2D, similitude3D, similitude,\n\
  affine2D, affine3D, affine, vectorfield2D, vectorfield3D, vectorfield, vector\n\
[-print] # print generated transformation\n\
[-no-print] # do not print generated transformation\n\
[-random] # random values\n\
[-identity] # identity transformation\n\
[-template %s] # reference image for geometry\n\
  the center of the image is a fixed point\n\
[-fixedpoint %lf %lf [%lf]] # fixed point\n";


static int _verbose_ = 1;


static void _ErrorParse( char *str, int flag );
static void _Parse( int argc, char *argv[], local_parameter *p );
static void _InitParam( local_parameter *par );
static char *_BaseName( char *p );


int main(int argc, char *argv[])
{
  local_parameter p;
  bal_transformation theTrsf;
  long int seedRandom = time(0);
 
  bal_image template;
 
  srandom( seedRandom );

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



  /***************************************************
   *
   * 
   *
   ***************************************************/
  BAL_InitTransformation( &theTrsf );
  theTrsf.type = p.transformation_type;



  switch ( p.transformation_value ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such value not handled yet for input transformation\n", program );
    return( -1 );

  case _RANDOM_ :
    
    switch ( p.transformation_type ) {
    default :
      if ( _verbose_ )
	fprintf( stderr, "%s: such type not handled yet for input transformation (identity case)\n", program );
      exit( -1 );

    case VECTORFIELD_2D :
    case VECTORFIELD_3D :
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

      if ( BAL_RandomTransformation( &theTrsf, p.transformation_type, (bal_image *)NULL ) != 1 ) {
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: unable to allocate transformation (matrix)\\n", program );
	exit( -1 );
      }
      break;

    }
    break;

  case _IDENTITY_ :

    switch ( p.transformation_type ) {
    default :
      if ( _verbose_ )
	fprintf( stderr, "%s: such type not handled yet for input transformation (identity case)\n", program );
      exit( -1 );

    case VECTORFIELD_2D :
    case VECTORFIELD_3D :
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

      if ( BAL_AllocTransformation( &theTrsf, p.transformation_type, (bal_image *)NULL ) != 1 ) {
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: unable to allocate transformation (matrix)\\n", program );
	exit( -1 );
      }
      break;
      
    }
    break;
    
  }
    

  switch ( p.transformation_value ) {

  case _IDENTITY_ :
    break;

  case _RANDOM_ :
  default :    
    
    /* transformation centering
     */
    if ( p.template_name != NULL ) {
      if ( BAL_ReadImage( &template, p.template_name, 0 ) != 1 ) {    
	fprintf( stderr, "%s: can not read '%s'\n", _BaseName(argv[0]), p.template_name );
	exit( -1 );
      }
      p.fixedpoint.x = ((template.ncols-1) * template.vx) / 2.0;
      p.fixedpoint.y = ((template.nrows-1) * template.vy) / 2.0;
      p.fixedpoint.z = ((template.nplanes-1) * template.vz) / 2.0;
      BAL_FreeImage( &template );
    }
    
    if ( p.fixedpoint.x > -9000 && p.fixedpoint.y > -9000 && p.fixedpoint.z > -9000 ) {
      
      switch ( p.transformation_type ) {
      case VECTORFIELD_2D :
      case VECTORFIELD_3D :
      default :
	if ( _verbose_ )
	  fprintf( stderr, "%s: such type not handled yet for input transformation (centering)\n", program );
	exit( -1 );
	
      case TRANSLATION_2D :
      case TRANSLATION_SCALING_2D :
      case RIGID_2D :
      case SIMILITUDE_2D :
      case AFFINE_2D :
	theTrsf.mat.m[3] = p.fixedpoint.x - theTrsf.mat.m[0] * p.fixedpoint.x - theTrsf.mat.m[1] * p.fixedpoint.y;
	theTrsf.mat.m[7] = p.fixedpoint.y - theTrsf.mat.m[4] * p.fixedpoint.x - theTrsf.mat.m[5] * p.fixedpoint.y;
	break;
      case TRANSLATION_3D :
      case TRANSLATION_SCALING_3D :
      case RIGID_3D :
      case SIMILITUDE_3D :
      case AFFINE_3D :
	theTrsf.mat.m[ 3] = p.fixedpoint.x - theTrsf.mat.m[0] * p.fixedpoint.x - theTrsf.mat.m[1] * p.fixedpoint.y - theTrsf.mat.m[ 2] * p.fixedpoint.z;
	theTrsf.mat.m[ 7] = p.fixedpoint.y - theTrsf.mat.m[4] * p.fixedpoint.x - theTrsf.mat.m[5] * p.fixedpoint.y - theTrsf.mat.m[ 6] * p.fixedpoint.z;
	theTrsf.mat.m[11] = p.fixedpoint.z - theTrsf.mat.m[8] * p.fixedpoint.x - theTrsf.mat.m[9] * p.fixedpoint.y - theTrsf.mat.m[10] * p.fixedpoint.z;
	break;
      }
      
    }

    break;
  }

  if ( p.print )
    BAL_PrintTransformation( stderr, &theTrsf, "generated transformation" );

  if ( p.restrsf_name != NULL ) {
    if ( BAL_WriteTransformation( &theTrsf, p.restrsf_name ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write '%s'\n", program, p.restrsf_name );
      BAL_FreeTransformation( &theTrsf );
      exit( -1 );
    }
  }
    


  BAL_FreeTransformation( &theTrsf );

  /*
  if (!createRandomTrsf(
              p.restrsf_name,
              p.template_name,
              p.fixedpoint,
              p.transformation_type,
              p.print))
  {
    fprintf( stderr, "%s: Failure.\n",program);
  }
*/

  exit( 0 );
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





static void _Parse( int argc, char *argv[], local_parameter *p )
{
  int i;
  int status;
  
  program = argv[0];
	
  for ( i=1; i<argc; i++ ) {
    

    if ( argv[i][0] == '-' ) {

      if ( strcmp ( argv[i], "-template") == 0
	   || (strcmp ( argv[i], "-t") == 0 && argv[i][2] == '\0')
	   || (strcmp ( argv[i], "-dims") == 0 && argv[i][5] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -reference", 0 );
	p->template_name = argv[i];
      }

      else if ( strcmp (argv[i], "-fixedpoint" ) == 0 ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -fixedpoint %lf", 0 );
	status = sscanf( argv[i], "%lf", &(p->fixedpoint.x) );
	if ( status <= 0 ) _ErrorParse( "parsing -fixedpoint %lf", 0 );
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -fixedpoint %lf %lf", 0 );
	status = sscanf( argv[i], "%lf", &(p->fixedpoint.y) );
	if ( status <= 0 ) _ErrorParse( "parsing -fixedpoint %lf %lf", 0 );
	i ++;
	if ( i >= argc) p->fixedpoint.z = 0;
	else {
	  status = sscanf( argv[i], "%lf", &(p->fixedpoint.z) );
	  if ( status <= 0 ) {
	    i--;
	    p->fixedpoint.z = 0;
	  }
	}  
      }
  
      /* transformation definition and computation 
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

      
      else if ( strcmp ( argv[i], "-random" ) == 0 
		|| (strcmp ( argv[i], "-rand" ) == 0 && argv[i][5] == '\0' ) )  {
	p->transformation_value = _RANDOM_;
      }
      else if ( strcmp ( argv[i], "-identity" ) == 0 
		|| (strcmp ( argv[i], "-id" ) == 0 && argv[i][3] == '\0' ) )  {
	p->transformation_value = _IDENTITY_;
      }


      else if ( strcmp ( argv[i], "-print" ) == 0 ) {
	p->print = 1;
      }
      else if ( strcmp ( argv[i], "-no-print" ) == 0 ) {
	p->print = 0;
      }

      else if ( (strcmp ( argv[i], "-verbose") == 0 && argv[i][8] == '\0')
		|| (strcmp ( argv[i], "-v") == 0 && argv[i][2] == '\0') ) {
	
      BAL_IncrementVerboseInBalImage(  );
      BAL_IncrementVerboseInBalTransformationTools(  );
      BAL_IncrementVerboseInBalTransformation(  );

    }

      else if ( strcmp ( argv[i], "--help" ) == 0 
		|| ( strcmp ( argv[i], "-help" ) == 0 && argv[i][5] == '\0' )
		|| ( strcmp ( argv[i], "--h" ) == 0 && argv[i][3] == '\0' )
		|| ( strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0' ) ) {
	_ErrorParse( NULL, 1 );
      }

      else {
	fprintf(stderr,"unknown option: '%s'\n",argv[i]);
      }

    }
    else {
      p->restrsf_name = argv[i];
    }
  }	
}





static void _InitParam( local_parameter *p ) 
{
  /* file names
     - images
     - transformations
  */
  p->restrsf_name = NULL;
  
  p->template_name = NULL;

  p->fixedpoint.x = -9999.0;
  p->fixedpoint.y = -9999.0;
  p->fixedpoint.z = -9999.0;

  p->transformation_type = AFFINE_3D;

  p->transformation_value = _RANDOM_;

  p->print = 0;
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


