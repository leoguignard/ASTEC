/*************************************************************************
 * blMaskDetectionResults.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2014, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Sam  8 fév 2014 19:21:38 CET
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
#include <bal-biolib-tools.h>






typedef struct local_parameter {

  /* file names
   */
  char *thedetection_name;
  char *resdetection_name;

  char *themask_name;

  /* format
   */
  char *thedetection_format;
  char *resdetection_format;
  int firstindex;
  int lastindex;

  char *themask_format;

} local_parameter;





/*------- Definition des fonctions statiques ----------*/
static void _ErrorParse( char *str, int flag );
static void _Parse( int argc, char *argv[], local_parameter *p );
static void _InitParam( local_parameter *par );




static char *program = NULL;

static char *usage = "%s %s [-mask |-m %s]\n\
 [-format %s] [-res-format %s] [-mask-format %s] -f[irst] %d -l[ast] %d\n\
 [-v] [-help]";

static char *detail = "\
-mask|-m %s # mask to be applied\n\
[-format %s] # format 'a la printf' of detection files to be processed\n\
                  # must contain one '%d'\n\
[-res-format %s]  # format 'a la printf' of result detection files\n\
                  # must contain one '%d'\n\
[-mask-format %s] # format 'a la printf' of mask files\n\
                  # must contain one '%d'\n\
[-first %d]        # first value of the index in the format\n\
[-last %d]         # last value of the index in the format\n\
 -v : mode verbose\n\
\n";






static int _process( char *thedetection_name,
		     char *resdetection_name,
		     char *themask_name );





int main( int argc, char *argv[] )
{
  local_parameter p;

  stringList theMaskFileList;
  stringList theDetectionFileList;
  stringList resDetectionFileList;
  char *theMaskName;
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
   


  /****************************************
   *
   * single file case
   *
   ****************************************/
  if ( p.thedetection_name != (char*)NULL
       && p.themask_name != (char*)NULL
       && p.resdetection_name != (char*)NULL ) {

    if ( _process( p.thedetection_name, p.resdetection_name,
		   p.themask_name ) != 0 ) {
      _ErrorParse( "unable to process single file case ...\n", 0);
    }
    
  }

  /****************************************
   *
   * format case
   *
   ****************************************/
  else  if ( p.thedetection_format != (char*)NULL
	     && ( p.themask_format != (char*)NULL
		  || p.themask_name != (char*)NULL )
	     && p.resdetection_format != (char*)NULL ) {

    if ( p.themask_format != (char*)NULL ) {
      initStringList( &theMaskFileList );
      
      if ( buildStringListFromFormat( p.themask_format, 
				      p.firstindex, p.lastindex, 
				      &theMaskFileList ) != 1 ) {
	_ErrorParse( "unable to build input transformation list from format\n", 0);
      }
    }
  
    initStringList( &theDetectionFileList );

    if ( buildStringListFromFormat( p.thedetection_format, 
				    p.firstindex, p.lastindex, 
				    &theDetectionFileList ) != 1 ) {
      if ( p.themask_format != (char*)NULL )
	freeStringList( &theMaskFileList );
      _ErrorParse( "unable to build input transformation list from format\n", 0);
    }
  
    initStringList( &resDetectionFileList );

    if ( buildStringListFromFormat( p.resdetection_format, 
				    p.firstindex, p.lastindex, 
				    &resDetectionFileList ) != 1 ) {
      freeStringList( &theDetectionFileList );
      if ( p.themask_format != (char*)NULL )
	freeStringList( &theMaskFileList );
      _ErrorParse( "unable to build input transformation list from format\n", 0);
    }
  
    theMaskName = p.themask_name;

    for ( i=0, n=p.firstindex; n<=p.lastindex; i++, n++ ) {
      if ( p.themask_format != (char*)NULL )
	theMaskName = theMaskFileList.data[i];
      if ( _process( theDetectionFileList.data[i],
		     resDetectionFileList.data[i],
		     theMaskName ) != 0 ) {
	if ( _verbose_ ) 
	  fprintf( stderr, " ... error when processing '%s' into '%s' with '%s'\n",
		   theDetectionFileList.data[i],
		   resDetectionFileList.data[i],
		   theMaskName );
	freeStringList( &resDetectionFileList );
	freeStringList( &theDetectionFileList );
	if ( p.themask_format != (char*)NULL )
	  freeStringList( &theMaskFileList );
	_ErrorParse( "unable to process format case ...\n", 0);
      }
    }


    freeStringList( &resDetectionFileList );
    freeStringList( &theDetectionFileList );
    if ( p.themask_format != (char*)NULL )
      freeStringList( &theMaskFileList );
  }

  /****************************************
   *
   * undefined case
   *
   ****************************************/
  else {
    _ErrorParse( "I don't know what to do ...\n", 0);
  }



  return( 1 );
}







static int _process( char *thedetection_name,
		     char *resdetection_name,
		     char *themask_name ) 
{
  char *proc = "_process";
  bal_blDetectionList thelist, reslist;
  bal_image theMask;
  int n;
  int x, y, z;
  
  if ( 0 ) {
    fprintf( stderr, "%s: processing '%s' into '%s' with '%s'\n",
	     proc, thedetection_name,
	     resdetection_name,
	     themask_name );
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
  if ( BAL_AllocBlDetectionList( &reslist, thelist.n ) != 0 ) {
    BAL_FreeBlDetectionList( &thelist );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate result detection\n", proc );
    return( -1 );
  }
  reslist.n = thelist.n;


  /* reading mask 
   */
  BAL_InitImage( &theMask, NULL, 0, 0, 0, 0, UCHAR );
  if ( BAL_ReadImage( &theMask, themask_name, 0 ) != 1 ) {
    _ErrorParse( "error when reading mask image\n", 0);
  }


  /* detection processing
   */
  switch ( theMask.type ) {
  default :
    BAL_FreeImage( &theMask );
    BAL_FreeBlDetectionList( &reslist );
    BAL_FreeBlDetectionList( &thelist );
    if ( _verbose_ )
      fprintf( stderr, "%s: mask image type not handled yet\n", 
	       proc );
    return( -1 );
    break;
  case UCHAR :
    {
      u8 ***theBuf = (u8***)theMask.array;
      for ( n=0, reslist.n=0; n<thelist.n; n++ ) {
	x = (int)(thelist.data[n].voxelcenter.x + 0.5);
	y = (int)(thelist.data[n].voxelcenter.y + 0.5);
	z = (int)(thelist.data[n].voxelcenter.z + 0.5);
	if ( x < 0 || x >= theMask.ncols ) continue;
	if ( y < 0 || y >= theMask.nrows ) continue;
	if ( z < 0 || z >= theMask.nplanes ) continue;
	if ( theBuf[z][y][x] == 0 ) continue;
	reslist.data[reslist.n] = thelist.data[n];
	reslist.n ++;
      }
    }
    break;
  }

  BAL_FreeImage( &theMask );
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

      /* mask related file names
       */

      if ( strcmp ( argv[i], "-mask") == 0
	   || (strcmp ( argv[i], "-m") == 0 && argv[i][2] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -mask", 0 );
	p->themask_name = argv[i];
      }
      
      else if ( strcmp ( argv[i], "-mask-format") == 0 ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -mask-format", 0 );
	p->themask_format = argv[i];
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

  p->themask_name = (char*)NULL;

  p->thedetection_format = (char*)NULL;
  p->resdetection_format = (char*)NULL;

  p->firstindex = 0;
  p->lastindex = 0;

  p->themask_format = (char*)NULL;

}
