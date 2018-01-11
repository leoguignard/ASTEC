/*************************************************************************
 * blApplyTrsfToDetectionResults.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2014, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Jeu 23 jan 2014 16:27:47 CET
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

#include <string-tools.h>

#include <bal-image.h>
#include <bal-transformation-tools.h>
#include <bal-biolib-tools.h>






typedef struct local_parameter {

  /* file names
   */
  char *thedetection_name;
  char *resdetection_name;

  char *therealtransformation_name;

  /* format
   */
  char *thedetection_format;
  char *resdetection_format;
  int firstindex;
  int lastindex;

  char *therealtransformation_format;

  /* template
  */
  char *theimage_name;
  char *resimage_name;

  int inversey;

} local_parameter;





/*------- Definition des fonctions statiques ----------*/
static void _ErrorParse( char *str, int flag );
static void _Parse( int argc, char *argv[], local_parameter *p );
static void _InitParam( local_parameter *par );




static char *program = NULL;

static char *usage = "%s %s [-transformation |-trsf %s]\n\
 [-template %s] [-res-template %s]\n\
 [-format %s] [-res-format %s] [-trsf-format %s] -f[irst] %d -l[ast] %d\n\
 [-inverse-y | -sx]\n\
 [-v] [-help]";

static char *detail = "\
-transformation|-trsf %s # transformation to be applied\n\
  # in 'real' coordinates.\n\
  # This transformation is with 'image' (ie applyTrsf) conventions\n\
  # and goes from 'result template' to 'template'.\n\
  # It will be inverted to be applied on points\n\
[-template %s] # template image corresponding to the input images\n\
  # (mandatory if transformation is given in real coordinates)\n\
[-res-template %s] # template image corresponding to the output images\n\
  # (mandatory if transformation is given in real coordinates)\n\
[-format %s] # format 'a la printf' of detection files to be processed\n\
                  # must contain one '%d'\n\
[-res-format %s]  # format 'a la printf' of result detection files\n\
                  # must contain one '%d'\n\
[-trsf-format %s] # format 'a la printf' of transformation files\n\
                  # must contain one '%d'\n\
[-first %d]        # first value of the index in the format\n\
[-last %d]         # last value of the index in the format\n\
[-inverse-y | -sx] # inverse l'axe des y  =  symetrie / x\n\
 -v : mode verbose\n\
\n";






static int _process( char *thedetection_name,
		     char *resdetection_name,
		     char *thetransformation_name,
		     bal_image *theTemplate,
		     bal_image *resTemplate );





int main( int argc, char *argv[] )
{
  local_parameter p;
  bal_image theTemplate;
  bal_image resTemplate;

  stringList theTransformationFileList;
  stringList theDetectionFileList;
  stringList resDetectionFileList;
  char *theTransformationName;
  int i, n;

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
  

  
  /* image templates
     (for real transformations)
  */
  if ( p.theimage_name == (char*)NULL ) {
    _ErrorParse( "image template(s) are mandatory for 'real' transformations\n", 0 );
  }
  BAL_InitImage( &theTemplate, NULL, 0, 0, 0, 0, UCHAR );
  BAL_InitImage( &resTemplate, NULL, 0, 0, 0, 0, UCHAR );
  if ( BAL_ReadImage( &theTemplate, p.theimage_name, 0 ) != 1 ) {
    _ErrorParse( "error when reading input template\n", 0);
  }
  if ( p.resimage_name == (char*)NULL ) {
    fprintf( stderr, "%s: uses input template as output template\n", program );
    resTemplate = theTemplate;
  }
  else {
    if ( BAL_ReadImage( &resTemplate, p.resimage_name, 0 ) != 1 ) {
      BAL_FreeImage( &theTemplate );
      _ErrorParse( "error when reading output template\n", 0);
    }
  }
  


  /****************************************
   *
   * single file case
   *
   ****************************************/
  if ( p.thedetection_name != (char*)NULL
       && p.therealtransformation_name != (char*)NULL
       && p.resdetection_name != (char*)NULL ) {

    if ( _process( p.thedetection_name, p.resdetection_name,
		   p.therealtransformation_name, &theTemplate, &resTemplate ) != 0 ) {
      if ( p.resimage_name != (char*)NULL ) BAL_FreeImage( &resTemplate );
      BAL_FreeImage( &theTemplate );
      _ErrorParse( "unable to process single file case ...\n", 0);
    }
    
  }

  /****************************************
   *
   * format case
   *
   ****************************************/
  else  if ( p.thedetection_format != (char*)NULL
	     && ( p.therealtransformation_format != (char*)NULL
		  || p.therealtransformation_name != (char*)NULL )
	     && p.resdetection_format != (char*)NULL ) {

    if ( p.therealtransformation_format != (char*)NULL ) {
      initStringList( &theTransformationFileList );
      
      if ( buildStringListFromFormat( p.therealtransformation_format, 
				      p.firstindex, p.lastindex, 
				      &theTransformationFileList ) != 1 ) {
	if ( p.resimage_name != (char*)NULL ) BAL_FreeImage( &resTemplate );
	BAL_FreeImage( &theTemplate );
	_ErrorParse( "unable to build input transformation file list from format\n", 0);
      }
    }
  
    initStringList( &theDetectionFileList );

    if ( buildStringListFromFormat( p.thedetection_format, 
				    p.firstindex, p.lastindex, 
				    &theDetectionFileList ) != 1 ) {
      if ( p.therealtransformation_format != (char*)NULL )
	freeStringList( &theTransformationFileList );
      if ( p.resimage_name != (char*)NULL ) BAL_FreeImage( &resTemplate );
      BAL_FreeImage( &theTemplate );
      _ErrorParse( "unable to build input detection file list from format\n", 0);
    }
  
    initStringList( &resDetectionFileList );

    if ( buildStringListFromFormat( p.resdetection_format, 
				    p.firstindex, p.lastindex, 
				    &resDetectionFileList ) != 1 ) {
      freeStringList( &theDetectionFileList );
      if ( p.therealtransformation_format != (char*)NULL )
	freeStringList( &theTransformationFileList );
      if ( p.resimage_name != (char*)NULL ) BAL_FreeImage( &resTemplate );
      BAL_FreeImage( &theTemplate );
      _ErrorParse( "unable to build output detection file  list from format\n", 0);
    }
  
    theTransformationName = p.therealtransformation_name;

    for ( i=0, n=p.firstindex; n<=p.lastindex; i++, n++ ) {
      if ( p.therealtransformation_format != (char*)NULL )
	theTransformationName = theTransformationFileList.data[i];
      if ( _process( theDetectionFileList.data[i],
		     resDetectionFileList.data[i],
		     theTransformationName,
		     &theTemplate, &resTemplate ) != 0 ) {
	if ( _verbose_ ) 
	  fprintf( stderr, " ... error when processing '%s' into '%s' with '%s'\n",
		   theDetectionFileList.data[i],
		   resDetectionFileList.data[i],
		   theTransformationName );
	freeStringList( &resDetectionFileList );
	freeStringList( &theDetectionFileList );
	if ( p.therealtransformation_format != (char*)NULL )
	  freeStringList( &theTransformationFileList );
	if ( p.resimage_name != (char*)NULL ) BAL_FreeImage( &resTemplate );
	BAL_FreeImage( &theTemplate );
	_ErrorParse( "unable to process format case ...\n", 0);
      }
    }


    freeStringList( &resDetectionFileList );
    freeStringList( &theDetectionFileList );
    if ( p.therealtransformation_format != (char*)NULL )
      freeStringList( &theTransformationFileList );
  }

  /****************************************
   *
   * undefined case
   *
   ****************************************/
  else {
    if ( p.resimage_name != (char*)NULL ) BAL_FreeImage( &resTemplate );
    BAL_FreeImage( &theTemplate );
    _ErrorParse( "I don't know what to do ...\n", 0);
  }



  if ( p.resimage_name != (char*)NULL ) BAL_FreeImage( &resTemplate );
  BAL_FreeImage( &theTemplate );

  return( 1 );
}







static int _process( char *thedetection_name,
		     char *resdetection_name,
		     char *thetransformation_name,
		     bal_image *theTemplate,
		     bal_image *resTemplate ) 
{
  char *proc = "_process";
  bal_blDetectionList thelist, reslist;
  bal_transformation tmpTrsf, invTrsf;
  int n;

  if ( 0 ) {
    fprintf( stderr, "%s: processing '%s' into '%s' with '%s'\n",
	     proc, thedetection_name,
	     resdetection_name,
	     thetransformation_name );
  }

  /* detection lists
   */
  BAL_InitBlDetectionList( &thelist );
  if ( BAL_ReadBlDetectionList( &thelist, thedetection_name ) != 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to read detection file '%s'\n", proc, thedetection_name );
    return( -1 );
  }
  
  BAL_InitBlDetectionList( &reslist );
  if ( BAL_AllocBlDetectionList( &reslist, thelist.n_allocated ) != 0 ) {
    BAL_FreeBlDetectionList( &thelist );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate result detection\n", proc );
    return( -1 );
  }
  reslist.n = thelist.n;



  /* transformations
   */
  BAL_InitTransformation( &tmpTrsf );
  BAL_InitTransformation( &invTrsf );
  

  /* read T_{image <- ref}
     ie T_{theTemplate <- resTemplate}
   */
  
  if ( BAL_ReadTransformation( &tmpTrsf, thetransformation_name ) != 1 ) {
    BAL_FreeBlDetectionList( &reslist );
    BAL_FreeBlDetectionList( &thelist );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to read '%s'\n", proc, thetransformation_name );
    return( -1 );
  }

  if ( BAL_AllocTransformation( &invTrsf, tmpTrsf.type, (bal_image*)NULL ) != 1 ) {
    BAL_FreeTransformation( &tmpTrsf );
    BAL_FreeBlDetectionList( &reslist );
    BAL_FreeBlDetectionList( &thelist );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate result transformation\n", proc );
    return( -1 );
  }

  /* calcul de T_{resTemplate <- theTemplate}
   */

  if ( BAL_InverseTransformation( &tmpTrsf, &invTrsf ) != 1 ) {
    BAL_FreeTransformation( &invTrsf );
    BAL_FreeTransformation( &tmpTrsf );
    BAL_FreeBlDetectionList( &reslist );
    BAL_FreeBlDetectionList( &thelist );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to invert transformation '%s'\n", proc, thetransformation_name );
    return( -1 );
  }
  
  if ( BAL_ChangeTransformationToVoxelUnit( resTemplate, theTemplate,
					    &invTrsf, &tmpTrsf ) != 1 ) {
    BAL_FreeTransformation( &invTrsf );
    BAL_FreeTransformation( &tmpTrsf );
    BAL_FreeBlDetectionList( &reslist );
    BAL_FreeBlDetectionList( &thelist );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to change units for transformation '%s'\n", 
	       proc, thetransformation_name );
    return( -1 );
  }


  BAL_FreeTransformation( &invTrsf );

  
  /* detection processing
   */
  for ( n=0; n<thelist.n; n++ ) {
    reslist.data[n] = thelist.data[n];
    if ( BAL_TransformPoint( &(thelist.data[n].voxelcenter), 
			     &(reslist.data[n].voxelcenter), &tmpTrsf ) != 1 ) {
      BAL_FreeTransformation( &tmpTrsf );
      BAL_FreeBlDetectionList( &reslist );
      BAL_FreeBlDetectionList( &thelist );
      if ( _verbose_ )
      fprintf( stderr, "%s: unable to transform center #%d of file '%s'\n", 
	       proc, n, thedetection_name );
    }
  }

  BAL_FreeTransformation( &tmpTrsf );
  BAL_FreeBlDetectionList( &thelist );
  


  /* result
   */
  
  if ( BAL_WriteBlDetectionList( &reslist, resdetection_name ) != 0 ) {
    BAL_FreeBlDetectionList( &reslist );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to write detection file '%s'\n", proc, resdetection_name );
    return( -1 );
  }
  
  BAL_FreeBlDetectionList( &reslist );

  return( 0 );
}










static void _Parse( int argc, char *argv[], local_parameter *p )
{
  int i, status;
  int inputisread = 0;
  int outputisread = 0;

  program = argv[0];
	
  for ( i=1; i<argc; i++ ) {
  
    if ( argv[i][0] == '-' ) {

      /* transformation related file names
       */

      if ( strcmp ( argv[i], "-transformation") == 0
	   || (strcmp ( argv[i], "-trsf") == 0 && argv[i][5] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -transformation", 0 );
	p->therealtransformation_name = argv[i];
      }
      
      else if ( strcmp ( argv[i], "-template") == 0
		|| (strcmp ( argv[i], "-t") == 0 && argv[i][2] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -template", 0 );
	p->theimage_name = argv[i];
      }
      else if ( strcmp ( argv[i], "-result-template") == 0
		|| strcmp ( argv[i], "-res-template") == 0
		|| (strcmp ( argv[i], "-res-t") == 0 && argv[i][6] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -res-template", 0 );
	p->resimage_name = argv[i];
      }

      else if ( strcmp ( argv[i], "-transformation-format") == 0
		|| strcmp ( argv[i], "-trsf-format") == 0 ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -transformation-format", 0 );
	p->therealtransformation_format = argv[i];
      }

      /* input related file names
       */
      else if ( strcmp ( argv[i], "-format") == 0 && argv[i][7] == '\0' ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -format", 0 );
	p->thedetection_format = argv[i];
      }


      /* output related file names
       */
      else if ( strcmp ( argv[i], "-result-format") == 0 
		||  strcmp ( argv[i], "-res-format") == 0  ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -result=format", 0 );
	p->resdetection_format = argv[i];
      }

      /* format related options
       */
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


      /* ...
       */

      else if ( (strcmp ( argv[i], "-inverse-y") == 0)
		|| (strcmp ( argv[i], "-sx") == 0 && argv[i][3] == '\0') ) {
	p->inversey = 1;
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
	p->thedetection_name = argv[i];
	inputisread = 1;
      }
      else if ( outputisread == 0 ) {
	p->resdetection_name = argv[i];
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
  p->thedetection_name = (char*)NULL;
  p->resdetection_name = (char*)NULL;

  p->therealtransformation_name = (char*)NULL;

  p->thedetection_format = (char*)NULL;
  p->resdetection_format = (char*)NULL;

  p->firstindex = 0;
  p->lastindex = 0;

  p->therealtransformation_format = (char*)NULL;

  p->theimage_name = (char*)NULL;
  p->resimage_name = (char*)NULL;

  p->inversey = 0;
}
