/*************************************************************************
 * createVectorTrsf.c -
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

#include <bal-createVectorTrsf.h>

typedef struct local_parameter {
  char *restrsf_name;

  char *template_image_name;

  bal_integerPoint dim;
  bal_doublePoint voxel;

  enumTypeTransfo transformation_type;

  enumVectorType vector_type;

  int print;

} local_parameter;


static char *program = NULL;

static char *usage = "%s\n\
 [-transformation-type|-transformation|-trsf-type %s] [-print | -no-print]\n\
 [-template %s] [-dim %d %d [%d]] [-voxel %lf %lf [%lf]]\n\
 [-vector-type %s]\n";

static char *detail = "\
[-transformation-type|-transformation|-trsf-type %s] # transformation type\n\
  vectorfield2D, vectorfield3D\n\
[-template %s]         # template image for geometry\n\
[-dim %d %d [%d]]      # template image dimensions\n\
[-voxel %lf %lf [%lf]] # voxel sizes of the template image\n\
[-vector-type %s]\n\
  sinus\n\
[-print] # print generated transformation\n\
[-no-print] # do not print generated transformation\n";


static int _verbose_ = 1;


static void _ErrorParse( char *str, int flag );
static void _Parse( int argc, char *argv[], local_parameter *p );
static void _InitParam( local_parameter *par );
static char *_BaseName( char *p );


int main(int argc, char *argv[])
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



  if (createVectorTrsf(
			p.restrsf_name,
			p.template_image_name,
			p.dim,
			p.voxel,
			p.transformation_type,
			p.vector_type,
			p.print,
			_verbose_
			))
    {
      fprintf( stderr, "%s: Failure.\n",program);
    }

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

      else if ( strcmp ( argv[i], "-vector-type" ) == 0 ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "-vector-type", 0 );
	if ( strcmp ( argv[i], "sinus2D" ) == 0 ) {
	  p->vector_type = SINUS_2D;
	}
	else if ( strcmp ( argv[i], "sinus3D" ) == 0 ) {
	  p->vector_type = SINUS_3D;
	}
	else if ( strcmp ( argv[i], "sinus" ) == 0 ) {
	  p->vector_type = SINUS_3D;
	}

	else {
	  fprintf( stderr, "unknown vector type: '%s'\n", argv[i] );
	  _ErrorParse( "-vector-type", 0 );
	}
      }

      else if ( strcmp ( argv[i], "-print" ) == 0 ) {
	p->print = 1;
      }
      else if ( strcmp ( argv[i], "-no-print" ) == 0 ) {
	p->print = 0;
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

  p->template_image_name = NULL;
  
  p->dim.x = 256;
  p->dim.y = 256;
  p->dim.z = 1;

  p->voxel.x = 1.0;
  p->voxel.y = 1.0;
  p->voxel.z = 1.0;

  p->transformation_type = VECTORFIELD_2D;

  p->vector_type = SINUS_2D;

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


